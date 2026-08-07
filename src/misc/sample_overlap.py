#!/usr/bin/env python3
"""
Sample overlap analysis across EPIC genetics studies.

Reads the Stage 3 final .psam files for all completed studies, builds a
per-study sample-ID set, then writes:

  1. An N x N overlap matrix TSV  — diagonal = study sample count,
     off-diagonal = shared samples between the two studies.
  2. A multi-study membership TSV — one row per sample that appears in
     more than one study, listing every study it belongs to.
  3. An UpSet plot PNG showing intersection sizes across studies.
  4. A human-readable summary printed to stdout.

Usage:
    python3 src/sample_overlap.py [--analysis-root DIR] [--out-dir DIR]
                                   [--no-plot]
                                   [--min-intersection INT]
                                   [--max-intersections INT]

The script auto-discovers .psam files at:
    <analysis-root>/<STUDY>/stage3/final/<STUDY>.psam

The UpSet plot requires matplotlib, upsetplot, and pandas.  If they are not
available the script runs normally and skips the plot with a warning.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent

# Stage 1 appends _R.. / _QC.. suffixes to replicate/QC sample IIDs.
SUFFIX_RE = re.compile(r"(_R[0-9]+|_QC[0-9]*)$")


class IdResolver:
    """Resolve a psam IID to its EPIC Idepic_Bio (participant) key.

    The genotype IID is a composite ``<Idepic>_<Idepic_Bio>`` for most samples and
    a single ``<Idepic_Bio>`` for others, so the same person can appear in two
    forms across studies. Resolving both to the Idepic_Bio unifies them. With a
    map (from the ``genetics_id.sas7bdat`` export) the resolution is authoritative;
    without one it parses the composite (second half = Idepic_Bio), which is
    correct for the standard fixed-width EPIC ID format. Either way, composite and
    single-field forms of the same participant are unified.
    """

    def __init__(self, id_map_path: Path | None):
        self.comp: dict[str, str] = {}
        self.by_bio: dict[str, str] = {}
        self.by_idp: dict[str, str] = {}
        if id_map_path is not None:
            with Path(id_map_path).open() as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                if not reader.fieldnames or "Idepic" not in reader.fieldnames or "Idepic_Bio" not in reader.fieldnames:
                    raise SystemExit(f"{id_map_path}: expected 'Idepic' and 'Idepic_Bio' columns")
                for row in reader:
                    a, b = row["Idepic"].strip(), row["Idepic_Bio"].strip()
                    self.comp[f"{a}_{b}"] = b
                    self.by_bio[b] = b
                    self.by_idp[a] = b

    def to_participant(self, iid: str) -> str:
        key = SUFFIX_RE.sub("", iid)
        for lookup in (self.comp, self.by_bio, self.by_idp):
            if key in lookup:
                return lookup[key]
        if len(key) == 29 and key[14] == "_":  # unmatched composite: take Idepic_Bio half
            return key[15:]
        return key


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_psam_ids(path: Path, resolver: IdResolver) -> set[str]:
    """Return the set of participant IDs (Idepic_Bio) from a PLINK2 .psam file.

    PLINK2 .psam format has a header line that begins with either:
      '#IID'  — no family ID column; the first column is the sample ID.
      '#FID'  — family ID is col 0, sample ID is col 1 (IID).

    Each IID is resolved to its authoritative Idepic_Bio so a participant is
    counted once across replicate suffixes and composite/single-field IID forms.
    """
    ids: set[str] = set()
    with path.open() as fh:
        header = fh.readline().rstrip("\n")
        cols = header.split("\t")
        if not cols:
            return ids
        if cols[0] == "#IID":
            iid_col = 0
        elif cols[0] == "#FID" and len(cols) > 1:
            iid_col = 1
        else:
            iid_col = 0
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) > iid_col:
                ids.add(resolver.to_participant(parts[iid_col]))
    return ids


def discover_psam_files(analysis_root: Path) -> dict[str, Path]:
    """Return {study_id: psam_path} for all Stage 3 final .psam files found.

    Stage 3 writes per-chromosome psam files (<STUDY>_chr<N>.psam); chr1 is used as
    the canonical representative since all chromosomes carry an identical sample list.
    """
    found: dict[str, Path] = {}
    for psam in sorted(analysis_root.glob("*/stage3/final/*_chr1.psam")):
        study = psam.parent.parent.parent.name
        found[study] = psam
    return found


# ---------------------------------------------------------------------------
# Matrix construction
# ---------------------------------------------------------------------------

def build_study_samples(psam_files: dict[str, Path], resolver: IdResolver) -> dict[str, set[str]]:
    study_samples: dict[str, set[str]] = {}
    for study, path in psam_files.items():
        ids = read_psam_ids(path, resolver)
        if ids:
            study_samples[study] = ids
        else:
            print(f"WARN: no samples read from {path}", file=sys.stderr)
    return study_samples


def overlap_matrix(
    studies: list[str],
    study_samples: dict[str, set[str]],
) -> list[list[int]]:
    n = len(studies)
    matrix: list[list[int]] = [[0] * n for _ in range(n)]
    for i, s1 in enumerate(studies):
        for j, s2 in enumerate(studies):
            matrix[i][j] = len(study_samples[s1] & study_samples[s2])
    return matrix


def multi_study_membership(
    study_samples: dict[str, set[str]],
) -> dict[str, list[str]]:
    """Return {iid: [study, ...]} for samples present in more than one study."""
    membership: dict[str, list[str]] = defaultdict(list)
    for study, ids in sorted(study_samples.items()):
        for iid in ids:
            membership[iid].append(study)
    return {iid: s for iid, s in membership.items() if len(s) > 1}


def compute_unique_count(
    study_samples: dict[str, set[str]],
) -> int:
    all_ids: set[str] = set()
    for ids in study_samples.values():
        all_ids |= ids
    return len(all_ids)


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_matrix_tsv(
    studies: list[str],
    matrix: list[list[int]],
    path: Path,
) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["Study"] + studies)
        for i, study in enumerate(studies):
            writer.writerow([study] + matrix[i])
    print(f"Overlap matrix written to:      {path}")


def write_membership_tsv(
    shared: dict[str, list[str]],
    path: Path,
) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["Idepic_Bio", "N_Studies", "Studies"])
        for participant, studies in sorted(shared.items()):
            writer.writerow([participant, len(studies), ",".join(studies)])
    print(f"Multi-study membership written to: {path}")


def print_summary(
    studies: list[str],
    study_samples: dict[str, set[str]],
    matrix: list[list[int]],
    shared: dict[str, list[str]],
) -> None:
    total_unique = compute_unique_count(study_samples)
    total_counted = sum(len(v) for v in study_samples.values())
    n_shared = len(shared)

    print()
    print("=" * 60)
    print(" EPIC Genetics — Sample Overlap Summary")
    print("=" * 60)
    print(f"  Studies analysed:          {len(studies)}")
    print(f"  Total unique samples:      {total_unique:,}")
    print(f"  Total sample-study pairs:  {total_counted:,}")
    print(f"  Samples in 2+ studies:     {n_shared:,}")
    print("=" * 60)
    print()

    print(f"  {'Study':<22}  {'N':>6}")
    print(f"  {'-'*22}  {'-'*6}")
    for study in studies:
        print(f"  {study:<22}  {len(study_samples[study]):>6,}")
    print()

    non_zero = [
        (studies[i], studies[j], matrix[i][j])
        for i in range(len(studies))
        for j in range(i + 1, len(studies))
        if matrix[i][j] > 0
    ]
    if non_zero:
        non_zero.sort(key=lambda x: -x[2])
        print(f"  Pairwise overlaps (non-zero only):")
        print(f"  {'Study A':<22}  {'Study B':<22}  {'Shared':>6}")
        print(f"  {'-'*22}  {'-'*22}  {'-'*6}")
        for s1, s2, count in non_zero:
            print(f"  {s1:<22}  {s2:<22}  {count:>6,}")
    else:
        print("  No pairwise sample overlap detected.")
    print()


# ---------------------------------------------------------------------------
# UpSet plot
# ---------------------------------------------------------------------------

def _build_combo_counts(
    studies: list[str],
    study_samples: dict[str, set[str]],
    shared: dict[str, list[str]],
) -> dict[frozenset, int]:
    """Return {frozenset(study_names): count} for every intersection."""
    combo_counts: dict[frozenset, int] = Counter()

    # Multi-study samples
    for studies_list in shared.values():
        combo_counts[frozenset(studies_list)] += 1

    # Single-study samples: total in study minus those also in other studies
    multi_in_study: Counter = Counter()
    for studies_list in shared.values():
        for s in studies_list:
            multi_in_study[s] += 1

    for study in studies:
        single_n = len(study_samples[study]) - multi_in_study[study]
        if single_n > 0:
            combo_counts[frozenset([study])] = single_n

    return combo_counts


def generate_upset_plot(
    studies: list[str],
    study_samples: dict[str, set[str]],
    shared: dict[str, list[str]],
    output_path: Path,
    min_intersection: int = 5,
    max_intersections: int = 25,
) -> None:
    """Generate and save an UpSet plot of sample overlap across studies."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mticker
        import pandas as pd
        from upsetplot import UpSet
    except ImportError as exc:
        print(
            f"WARN: UpSet plot skipped — missing package ({exc}). "
            "Install with: pip install upsetplot matplotlib pandas",
            file=sys.stderr,
        )
        return

    combo_counts = _build_combo_counts(studies, study_samples, shared)
    if not combo_counts:
        print("WARN: No intersections to plot — plot skipped.", file=sys.stderr)
        return

    # Sort studies by total N descending so the largest studies appear at top
    study_order = sorted(studies, key=lambda s: -len(study_samples[s]))

    # Build the data Series from EVERY intersection — the full, disjoint
    # partition of all samples. This is essential: upsetplot derives the
    # per-study set-size totals from the complete input, so passing the full
    # partition makes the totals bar equal the true study N. (The previous
    # version pre-filtered the data before handing it to upsetplot, which made
    # the totals undercount — e.g. Brea_01_Erneg showed 930 instead of 987.)
    tuples = [tuple(s in combo for s in study_order) for combo in combo_counts]
    counts = list(combo_counts.values())

    mi = pd.MultiIndex.from_tuples(tuples, names=study_order)
    data = pd.Series(counts, index=mi)

    # Remove all-False rows (shouldn't exist, but guard against it)
    data = data[data.index.to_frame().any(axis=1)]

    # Limit only the *displayed* intersection bars via `min_subset_size` so the
    # chart stays readable; totals are unaffected because they come from the
    # full data above. Raise the effective minimum to the size at rank
    # `max_intersections` when there are more qualifying intersections than that.
    sizes_desc = sorted(combo_counts.values(), reverse=True)
    effective_min = min_intersection
    if len(sizes_desc) > max_intersections:
        effective_min = max(effective_min, sizes_desc[max_intersections - 1])

    # Never drop a study entirely: every study must keep at least its singleton
    # (unique-sample) column. Each study's unique count is the size of its
    # single-study combo; cap the threshold at the smallest of these so the
    # `max_intersections` rule can't push it above a small study's singleton.
    singleton_sizes = [
        combo_counts[frozenset([s])]
        for s in study_order
        if frozenset([s]) in combo_counts
    ]
    if singleton_sizes:
        effective_min = min(effective_min, min(singleton_sizes))

    # --- Figure ---
    # Size the canvas to the content: width scales with the number of displayed
    # intersection columns, height with the number of study rows. A fixed width
    # clips the matrix/bars once there are many columns, so compute it instead.
    n_cols = sum(1 for v in combo_counts.values() if v >= effective_min)
    fig_w = max(36.0, 1.6 * n_cols + 14.0)
    fig_h = max(12.0, 0.55 * len(study_order) + 5.0)
    fig = plt.figure(figsize=(fig_w, fig_h))
    upset = UpSet(
        data,
        subset_size="sum",
        show_counts="%d",
        sort_by="cardinality",
        min_subset_size=effective_min,
        totals_plot_elements=3,
    )
    upset.style_subsets(
        present=list(study_order),
        facecolor="#3b82f6",
    )
    axes = upset.plot(fig)

    # Style the totals (set-size) bar chart
    ax_totals = axes["totals"]
    ax_totals.set_xlabel("Samples per study", fontsize=9)
    ax_totals.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

    # Style the intersection bar chart
    ax_inter = axes["intersections"]
    ax_inter.set_ylabel("Intersection size", fontsize=9)
    ax_inter.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

    # Caption (bottom-left): note that the figure is simplified — intersections
    # smaller than the display threshold are not drawn.
    fig.text(
        0.01, 0.01,
        f"simplified figure; intersections <{int(effective_min)} excluded",
        ha="left", va="bottom", fontsize=9, style="italic", color="#666666",
    )

    # Choose a DPI that keeps the rendered image within matplotlib's hard
    # 65,536-px-per-dimension limit, so requesting "show every intersection"
    # (a very wide canvas) cannot crash the save.
    dpi = 150
    max_px = 60000
    longest_in = max(fig_w, fig_h)
    if longest_in * dpi > max_px:
        dpi = max(50, int(max_px / longest_in))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"UpSet plot saved to:            {output_path}  ({fig_w:.0f}x{fig_h:.0f} in @ {dpi} dpi)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute sample overlap across EPIC genetics Stage 3 studies."
    )
    parser.add_argument(
        "--analysis-root",
        default=str(REPO_ROOT / "analysis"),
        help="Path to the analysis root directory (default: <repo>/analysis)",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help=(
            "Directory for output files (default: <analysis-root>/report). "
            "Created if it does not exist."
        ),
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip UpSet plot generation.",
    )
    parser.add_argument(
        "--min-intersection",
        type=int,
        default=5,
        help="Minimum intersection size to display in the UpSet plot (default: 5).",
    )
    parser.add_argument(
        "--max-intersections",
        type=int,
        default=25,
        help="Maximum number of multi-study intersections to show (default: 25).",
    )
    parser.add_argument(
        "--id-map",
        default=None,
        help=(
            "TSV with Idepic and Idepic_Bio columns (exported by 003-data-epic.R). "
            "Resolves each IID to its authoritative Idepic_Bio so composite and "
            "single-field forms of the same person are unified. Without it, counts "
            "over-count unique participants and under-count overlaps."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    analysis_root = Path(args.analysis_root).resolve()

    if not analysis_root.is_dir():
        raise SystemExit(f"ERROR: analysis root not found: {analysis_root}")

    out_dir = Path(args.out_dir).resolve() if args.out_dir else analysis_root / "report"
    out_dir.mkdir(parents=True, exist_ok=True)

    psam_files = discover_psam_files(analysis_root)
    if not psam_files:
        raise SystemExit(
            f"ERROR: no Stage 3 final .psam files found under {analysis_root}"
        )
    print(f"Found {len(psam_files)} studies with Stage 3 final .psam files.")

    if args.id_map:
        resolver = IdResolver(Path(args.id_map).resolve())
    else:
        resolver = IdResolver(None)
        print(
            "NOTE: no --id-map given; resolving Idepic_Bio by parsing the composite "
            "IID (correct for the standard fixed-width EPIC ID format). Pass "
            "--id-map for an authoritative, source-validated resolution.",
            file=sys.stderr,
        )

    study_samples = build_study_samples(psam_files, resolver)
    studies = sorted(study_samples.keys(), key=str.lower)

    matrix = overlap_matrix(studies, study_samples)
    shared = multi_study_membership(study_samples)

    write_matrix_tsv(studies, matrix, out_dir / "sample_overlap_matrix.tsv")
    if shared:
        write_membership_tsv(shared, out_dir / "sample_overlap_membership.tsv")
    else:
        print("No samples shared across studies — membership file not written.")

    if not args.no_plot:
        generate_upset_plot(
            studies,
            study_samples,
            shared,
            out_dir / "sample_overlap_upset.png",
            min_intersection=args.min_intersection,
            max_intersections=args.max_intersections,
        )

    print_summary(studies, study_samples, matrix, shared)

    total_unique = compute_unique_count(study_samples)
    summary = {
        "n_studies": len(studies),
        "total_unique": total_unique,
        "total_pairs": sum(len(v) for v in study_samples.values()),
        "n_shared": len(shared),
    }
    summary_json_path = out_dir / "sample_overlap_summary.json"
    summary_json_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Summary JSON written to {summary_json_path}")


if __name__ == "__main__":
    main()
