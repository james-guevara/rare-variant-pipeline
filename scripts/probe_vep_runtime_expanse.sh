#!/usr/bin/env bash
# Slurm diagnostic for the packaged Perl/module paths in Expanse's VEP image.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --partition=ind-shared

set -euo pipefail
module load singularitypro/4.1.2
image=${VEP_IMAGE:?VEP_IMAGE is required}
singularity exec "$image" sh -c '
  command -v perl
  command -v vep
  find /opt /usr/local -path "*/Bio/EnsEMBL/Attribute.pm" 2>/dev/null | head
'
