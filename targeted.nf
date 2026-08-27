#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.manifest = null
params.bindings = null
params.targeted_container = null
params.preflight_only = false
params.run_root = null

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

workflow {
    if (!params.manifest || !params.bindings || !params.targeted_container) {
        error 'Required: --manifest, --bindings, and --targeted_container'
    }
    manifests = Channel.of(tuple(file(params.manifest), file(params.bindings)))
    TARGETED_CHROMOSOME(manifests)
}
