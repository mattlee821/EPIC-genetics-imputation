process SAMPLE_MANIFEST {
    tag "sample_manifest"
    publishDir "${params.dest_root}", mode: 'copy', overwrite: true

    input:
    path psams

    output:
    path("sample-manifest.tsv")

    script:
    """
    \$PYTHON3_BIN "${params.project_root}/pipeline_stage4/bin/build_sample_manifest.py" \\
      --id-map "${params.idepic_map}" \\
      --out sample-manifest.tsv \\
      ${psams}
    """
}
