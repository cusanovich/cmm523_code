#!/bin/bash
# ============================================================================
# CMM 523 -- render a tutorial as a batch job
#
#   sbatch render_tutorial.sbatch 09_wnn.qmd
#
# For the tutorials that take hours. The Signac, RNA+ATAC and WNN tutorials run
# with caching disabled (knitr's cache cannot hold objects over 2 GB), so every
# render is a full run.
#
# Watch:  squeue --me   /   tail -f render_<jobid>.log
# ============================================================================

#SBATCH --job-name=cmm523_render
#SBATCH --account=cusanovichlab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=10:00:00
#SBATCH --output=render_%j.log

set -e

QMD=$1
if [ -z "${QMD}" ]; then
    echo "Usage: sbatch render_tutorial.sbatch <tutorial.qmd>"
    exit 1
fi

SIF=/groups/darrenc/cmm_523/containers/cmm523_full.sif
export CMM523_NETID=darrenc
export CMM523_CACHE=/xdisk/darrenc/darrenc/cmm523_cache

cd "${SLURM_SUBMIT_DIR}"

echo "=========================================="
echo "Node:     $(hostname)"
echo "Started:  $(date)"
echo "Tutorial: ${QMD}"
echo "=========================================="

if [ ! -f "${QMD}" ]; then
    echo "ERROR: ${QMD} not found in $(pwd)"
    exit 1
fi

apptainer exec -B /xdisk -B /groups "${SIF}" quarto render "${QMD}"

echo "=========================================="
echo "Finished: $(date)"
ls -la "${QMD%.qmd}.md"
echo "=========================================="
