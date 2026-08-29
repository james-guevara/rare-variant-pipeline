#!/usr/bin/env nextflow
process TARGETED_CHROMOSOME {
    tag { "${science_manifest.simpleName}" }
    container params.targeted_container
    cpus 4
    memory '16 GB'
    time '2h'

    input:
    tuple path(science_manifest), path(environment_bindings)

    output:
    path 'execution.receipt.json'
    path 'targeted-output'

    script:
    def preflight = params.preflight_only ? '--preflight-only' : ''
    def runRoot = params.run_root ? "--run-root ${params.run_root}" : ''
    """
    set -euo pipefail
    python /opt/rvp/scripts/run_targeted_manifest.py \
      --manifest ${science_manifest} \
      --bindings ${environment_bindings} \
      ${runRoot} \
      ${preflight}
    python /opt/rvp/scripts/package_targeted_outputs.py \
      --science-manifest ${science_manifest} \
      --bindings ${environment_bindings} \
      ${runRoot} \
      --output-dir targeted-output \
      ${preflight}
    python - <<'PY' > execution.receipt.json
    import json
    print(json.dumps({
        "manifest": "${science_manifest.name}",
        "bindings": "${environment_bindings.name}",
        "preflight_only": "${params.preflight_only}" == "true",
        "status": "SUCCEEDED",
    }, sort_keys=True))
    PY
    """
}

workflow TARGETED_MANIFEST_WORKFLOW {
    take:
    manifest_bindings

    main:
    TARGETED_CHROMOSOME(manifest_bindings)

    emit:
    execution_receipt = TARGETED_CHROMOSOME.out[0]
    chromosome_outputs = TARGETED_CHROMOSOME.out[1]
}
