# Source from Nextflow beforeScript, after loading Singularity.
# Relocate this exact immutable image only; preserve Nextflow's image identity.
set -euo pipefail
source_image=${1:?Expected source SIF path}
expected=${2:?Expected SHA-256}
[[ "$source_image" == /* && -f "$source_image" && "$expected" =~ ^[0-9a-f]{64}$ ]] || return 1
runtime_root=/scratch/$USER/job_${SLURM_JOB_ID:?}/rare-container
mkdir -p "$runtime_root/bin"
real_singularity=$(command -v singularity)
cp "$source_image" "$runtime_root/image.sif.partial"
printf '%s  %s\n' "$expected" "$runtime_root/image.sif.partial" | sha256sum -c -
mv "$runtime_root/image.sif.partial" "$runtime_root/image.sif"
# Embed paths because Nextflow deliberately clears most environment variables.
python3 - "$runtime_root" "$source_image" "$real_singularity" <<'PY'
import sys,shlex,pathlib
root,source,binary=sys.argv[1:]
p=pathlib.Path(root)/'bin/singularity'
p.write_text('#!/bin/bash\nset -euo pipefail\nargs=()\nfor arg in "$@"; do\n if [[ "$arg" == '+shlex.quote(source)+' ]]; then arg='+shlex.quote(root+'/image.sif')+'; fi\n args+=("$arg")\ndone\nexec '+shlex.quote(binary)+' "${args[@]}"\n')
p.chmod(0o700)
PY
export PATH="$runtime_root/bin:$PATH"
singularity exec "$source_image" python -c 'import gzip,zlib; print("Verified node-local container: gzip/zlib imports passed")'
printf 'source=%s\nsha256=%s\nlocal=%s\nnode=%s\n' "$source_image" "$expected" "$runtime_root/image.sif" "$HOSTNAME" > container-localization.txt
