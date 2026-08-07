#!/usr/bin/env python3
"""Build a combined sample manifest (study, sample ID, EPIC IDs, sex, phenotype,
overlap flag) from the finalised per-study PLINK2 .psam files, for sharing with
collaborators to check which of their samples overlap the genetics data.

Called by the Stage 4 SAMPLE_MANIFEST process with the canonical per-study
.psam files staged into the working directory. One row per (study, sample).
Study names are taken from the .psam filename (``<study>_chr1.psam``). Column
positions are read from each file's header, so it is robust to which optional
columns PLINK2 wrote.

The genotype sample IID is a composite ``<Idepic>_<Idepic_Bio>`` for most samples
and a single ``<Idepic_Bio>`` for others. Rather than string-parse it, the two
EPIC IDs are looked up from the authoritative id map (``--id-map``: a TSV with
``Idepic`` and ``Idepic_Bio`` columns, exported from genetics_id.sas7bdat by
003-data-epic.R). This resolves both IDs correctly even for single-field samples.

Columns:
  study        study ID
  IID          sample ID exactly as stored in the .psam (the plink2 --keep key)
  Idepic       EPIC participant (person) ID, from the id map
  Idepic_Bio   EPIC biosample ID, from the id map — the key for matching cohorts
  sex          PLINK sex code as stored (1=male, 2=female, 0/NA=unknown)
  pheno        phenotype as stored (PLINK 1=control, 2=case, NA=missing)
  overlap_keep participants in >1 study are kept once, in the study with the
               smallest total N (ties broken by study name): TRUE in that study,
               FALSE in the others; NA for participants in a single study.
               Keyed on Idepic_Bio, so the same person in different IID forms
               across studies is one participant.

Usage:
  build_sample_manifest.py --id-map EPIC_Idepic_map.tsv --out sample-manifest.tsv \\
      <study>_chr1.psam [...]

Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

SUFFIX_RE = re.compile(r"(_R[0-9]+|_QC[0-9]*)$")
STUDY_FROM_PSAM_RE = re.compile(r"_chr[0-9XY]+\.psam$", re.IGNORECASE)


def study_from_psam(path: Path) -> str:
    name = STUDY_FROM_PSAM_RE.sub("", path.name)
    return name[:-5] if name.endswith(".psam") else name  # fallback: strip .psam


def parse_psam(path: Path) -> tuple[list[str], int, int, int | None]:
    """Return (data_lines, iid_idx, sex_idx, pheno_idx) for a PLINK2 .psam."""
    lines = path.read_text().splitlines()
    header = next((ln for ln in lines if ln.startswith("#")), None)
    if header is None:
        raise ValueError(f"{path}: no '#' header line (not a PLINK2 .psam?)")
    cols = header.lstrip("#").split()
    idx = {name: i for i, name in enumerate(cols)}

    if "IID" not in idx:
        raise ValueError(f"{path}: header has no IID column ({cols})")
    if "SEX" not in idx:
        raise ValueError(f"{path}: header has no SEX column ({cols})")

    reserved = {"FID", "IID", "PAT", "MAT", "SEX"}
    pheno_idx = idx.get("PHENO1", idx.get("PHENO"))
    if pheno_idx is None:
        pheno_idx = next((i for i, name in enumerate(cols) if name not in reserved), None)

    data = [ln for ln in lines if ln and not ln.startswith("#")]
    return data, idx["IID"], idx["SEX"], pheno_idx


class IdResolver:
    """Resolve a psam IID to authoritative (Idepic, Idepic_Bio) via the id map."""

    def __init__(self, id_map_path: Path):
        self.by_composite: dict[str, tuple[str, str]] = {}
        self.by_bio: dict[str, tuple[str, str]] = {}
        self.by_idepic: dict[str, tuple[str, str]] = {}
        with id_map_path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if not reader.fieldnames or "Idepic" not in reader.fieldnames or "Idepic_Bio" not in reader.fieldnames:
                raise ValueError(f"{id_map_path}: expected 'Idepic' and 'Idepic_Bio' columns, got {reader.fieldnames}")
            for row in reader:
                a, b = row["Idepic"].strip(), row["Idepic_Bio"].strip()
                self.by_composite[f"{a}_{b}"] = (a, b)
                self.by_bio[b] = (a, b)
                self.by_idepic[a] = (a, b)

    def resolve(self, iid: str) -> tuple[str, str]:
        key = SUFFIX_RE.sub("", iid)
        for lookup in (self.by_composite, self.by_bio, self.by_idepic):
            if key in lookup:
                return lookup[key]
        # Not in the id map: best-effort parse of the composite, else Idepic unknown.
        if len(key) == 29 and key[14] == "_":
            return key[:14], key[15:]
        return "NA", key


def overlap_keep_flags(rows, per_study_n):
    """Map (study, iid) -> 'TRUE'/'FALSE'/'NA' for the overlap_keep column.

    Participants (keyed on Idepic_Bio) in more than one study are kept once, in
    the study with the smallest total N (ties broken by study name). NA for a
    participant in a single study.
    """
    part_studies = defaultdict(set)
    for study, _iid, _idepic, idepic_bio, _sex, _pheno in rows:
        part_studies[idepic_bio].add(study)

    keep_study = {}  # Idepic_Bio -> chosen study (only for multi-study participants)
    for idepic_bio, studies in part_studies.items():
        if len(studies) > 1:
            keep_study[idepic_bio] = min(studies, key=lambda s: (per_study_n[s], s))

    flags = {}
    for study, iid, _idepic, idepic_bio, _sex, _pheno in rows:
        if idepic_bio not in keep_study:
            flags[(study, iid)] = "NA"
        else:
            flags[(study, iid)] = "TRUE" if study == keep_study[idepic_bio] else "FALSE"
    return flags


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("psams", nargs="+", type=Path, help="Per-study PLINK2 .psam files")
    ap.add_argument("--id-map", type=Path, required=True,
                    help="TSV with Idepic and Idepic_Bio columns (from 003-data-epic.R)")
    ap.add_argument("--out", type=Path, required=True, help="Output TSV path")
    args = ap.parse_args()

    resolver = IdResolver(args.id_map)

    rows: list[tuple[str, str, str, str, str, str]] = []
    per_study_n: Counter[str] = Counter()
    n_unresolved = 0
    for psam in args.psams:
        study = study_from_psam(psam)
        data, iid_i, sex_i, pheno_i = parse_psam(psam)
        for line in data:
            fields = line.split()
            iid = fields[iid_i]
            idepic, idepic_bio = resolver.resolve(iid)
            n_unresolved += idepic == "NA"
            sex = fields[sex_i] if sex_i < len(fields) else ""
            pheno = fields[pheno_i] if pheno_i is not None and pheno_i < len(fields) else ""
            rows.append((study, iid, idepic, idepic_bio, sex, pheno))
        per_study_n[study] = len(data)

    flags = overlap_keep_flags(rows, per_study_n)

    rows.sort(key=lambda r: (r[0], r[1]))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["study", "IID", "Idepic", "Idepic_Bio", "sex", "pheno", "overlap_keep"])
        for study, iid, idepic, idepic_bio, sex, pheno in rows:
            writer.writerow([study, iid, idepic, idepic_bio, sex, pheno, flags[(study, iid)]])

    print(f"Wrote {len(rows)} samples across {len(per_study_n)} studies -> {args.out}", file=sys.stderr)
    if n_unresolved:
        print(f"WARNING: {n_unresolved} sample(s) not found in id map (Idepic=NA)", file=sys.stderr)
    for study, n in sorted(per_study_n.items()):
        print(f"  {study:16s} {n:>7d}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
