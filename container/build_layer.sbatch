#!/bin/bash
# ============================================================================
# CMM 523 -- build any container layer as a batch job
#
#   sbatch build_layer.sbatch cmm523_base.def  cmm523_base.sif   # ~90 min
#   sbatch build_layer.sbatch cmm523_tools.def cmm523_2026.sif   # ~15 min
#
# Watch:  squeue --me   /   tail -f build_<jobid>.log
# ============================================================================

#SBATCH --job-name=cmm523_build
#SBATCH --account=cusanovichlab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=build_%j.log

set -e

DEF_NAME=$1
SIF_NAME=$2

if [ -z "${DEF_NAME}" ] || [ -z "${SIF_NAME}" ]; then
    echo "Usage: sbatch build_layer.sbatch <recipe.def> <output.sif>"
    exit 1
fi

DEF=${SLURM_SUBMIT_DIR}/${DEF_NAME}
SIF=/xdisk/darrenc/cmm_523/containers/${SIF_NAME}

export APPTAINER_CACHEDIR=/xdisk/darrenc/cmm_523/.apptainer_cache
export APPTAINER_TMPDIR=/xdisk/darrenc/cmm_523/.apptainer_tmp
mkdir -p "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}" "$(dirname ${SIF})"

echo "=========================================="
echo "Node:    $(hostname)"
echo "Started: $(date)"
echo "Recipe:  ${DEF}"
echo "Output:  ${SIF}"
echo "=========================================="

if [ ! -f "${DEF}" ]; then
    echo "ERROR: recipe not found: ${DEF}"
    exit 1
fi

if ! command -v apptainer > /dev/null 2>&1; then
    echo "ERROR: apptainer not found on this node."
    echo "PATH is: ${PATH}"
    exit 1
fi

apptainer --version

# Move any previous image aside rather than clobbering a working one.
if [ -f "${SIF}" ]; then
    mv "${SIF}" "${SIF}.$(date +%Y%m%d_%H%M%S).bak"
    echo "Moved previous image aside."
fi

apptainer build "${SIF}" "${DEF}"

echo "=========================================="
echo "Finished: $(date)"
ls -lh "${SIF}"
echo "=========================================="

apptainer exec -B /xdisk "${SIF}" Rscript -e '
  cat(R.version.string, "\n")
  cat("Packages installed:", nrow(installed.packages()), "\n")
'
