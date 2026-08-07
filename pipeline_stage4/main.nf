#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MASTER_REPORT   } from './modules/report.nf'
include { FINALISE_STUDY  } from './modules/finalise.nf'
include { SAMPLE_MANIFEST } from './modules/manifest.nf'

workflow {

    def included_studies = params.study == 'all' ? [] : params.study.split(',').collect { it.trim() }

    // ── Per-chromosome pgen/pvar files published by stage 3 FINALIZE_CHROM ──────
    // Path: ${stage3_root}/${study}/stage3/final/${study}_chr*.pgen
    ch_pgens = Channel.fromPath("${params.stage3_root}/*/stage3/final/*.pgen")
        .map { f -> tuple(f.parent.parent.parent.name, f) }
        .filter { study, f -> params.study == 'all' || included_studies.contains(study) }
        .groupTuple()

    ch_pvars = Channel.fromPath("${params.stage3_root}/*/stage3/final/*.pvar")
        .map { f -> tuple(f.parent.parent.parent.name, f) }
        .filter { study, f -> params.study == 'all' || included_studies.contains(study) }
        .groupTuple()

    // All per-chr psams are identical; chr1 is the canonical representative.
    // One branch feeds FINALISE_STUDY (per study), the other the combined manifest.
    Channel.fromPath("${params.stage3_root}/*/stage3/final/*_chr1.psam")
        .map { f -> tuple(f.parent.parent.parent.name, f) }
        .filter { study, f -> params.study == 'all' || included_studies.contains(study) }
        .multiMap { study, psam ->
            for_finalise: tuple(study, psam)
            for_manifest: psam
        }
        .set { ch_psam_split }

    // ── QC outputs from stage 3 SAMPLE_REVIEW_SUMMARY and HET_PCA_QC/KING_QC ───
    // Exclude files: ${study}/stage3/report/flags/
    ch_related = Channel.fromPath("${params.stage3_root}/*/stage3/report/flags/*.related.exclude")
        .map { f -> tuple(f.parent.parent.parent.parent.name, f) }
        .filter { study, f -> params.study == 'all' || included_studies.contains(study) }

    ch_ancestry = Channel.fromPath("${params.stage3_root}/*/stage3/report/flags/*.ancestry.exclude")
        .map { f -> tuple(f.parent.parent.parent.parent.name, f) }
        .filter { study, f -> params.study == 'all' || included_studies.contains(study) }

    ch_hwe = Channel.fromPath("${params.stage3_root}/*/stage3/report/flags/*.hwe.exclude")
        .map { f -> tuple(f.parent.parent.parent.parent.name, f) }
        .filter { study, f -> params.study == 'all' || included_studies.contains(study) }

    // Review files: ${study}/stage3/sample_review/
    ch_het = Channel.fromPath("${params.stage3_root}/*/stage3/sample_review/*_het.het")
        .map { f -> tuple(f.parent.parent.parent.name, f) }
        .filter { study, f -> params.study == 'all' || included_studies.contains(study) }

    ch_kin0 = Channel.fromPath("${params.stage3_root}/*/stage3/sample_review/*_king.kin0")
        .map { f -> tuple(f.parent.parent.parent.name, f) }
        .filter { study, f -> params.study == 'all' || included_studies.contains(study) }

    ch_eigenvec = Channel.fromPath("${params.stage3_root}/*/stage3/sample_review/*_pca.eigenvec")
        .map { f -> tuple(f.parent.parent.parent.name, f) }
        .filter { study, f -> params.study == 'all' || included_studies.contains(study) }

    ch_eigenval = Channel.fromPath("${params.stage3_root}/*/stage3/sample_review/*_pca.eigenval")
        .map { f -> tuple(f.parent.parent.parent.name, f) }
        .filter { study, f -> params.study == 'all' || included_studies.contains(study) }

    // ── Per-study stage HTML reports ─────────────────────────────────────────────
    // stage3 HTML is written by run_stage3_reports.py (called at end of 006_stage3.sh)
    ch_stage3_html = Channel.fromPath("${params.stage3_root}/*/stage3/report/report-stage3.html")
        .map { f -> tuple(f.parent.parent.parent.name, f) }
        .filter { study, f -> params.study == 'all' || included_studies.contains(study) }
        .multiMap { study, html ->
            for_master:   tuple(study, html)
            for_finalise: tuple(study, html)
        }
        .set { ch_stage3_split }

    // ── MASTER_REPORT: generate per-study cross-stage master HTML ────────────────
    // Uses stage3_html as cache sentinel; reads full report tree from analysis_root
    // internally. stage2_html is looked up from analysis_root in FINALISE_STUDY.
    ch_master = MASTER_REPORT(
        ch_stage3_split.for_master
    )

    // ── FINALISE_STUDY: build deliverable archive + copy HTML + review files ─────
    ch_finalise_input = ch_pgens
        .join(ch_pvars)
        .join(ch_psam_split.for_finalise)
        .join(ch_related)
        .join(ch_ancestry)
        .join(ch_hwe)
        .join(ch_stage3_split.for_finalise)
        .join(ch_master)
        .join(ch_het)
        .join(ch_kin0)
        .join(ch_eigenvec)
        .join(ch_eigenval)

    FINALISE_STUDY(ch_finalise_input)

    // ── SAMPLE_MANIFEST: combined study/IID/sex/pheno table for collaborators ────
    SAMPLE_MANIFEST(ch_psam_split.for_manifest.collect())
}
