#!/bin/bash
# ============================================================================
# CMM 523 -- render every tutorial, in order
#
#   cd ~/ccsgab_projects/cmm523_code/tutorials
#   sbatch render_all.sbatch
#
# Order matters: tutorial 9 reads the object that tutorial 8 saves, so 8 has to
# run first. The rest are independent.
#
# Roughly 5-6 hours in total if the caches are empty. Tutorials 7, 8 and 9 do
# not cache at all (knitr's cache cannot hold an object over 2 GB), so those
# three always run in full.
#
# Watch:  squeue --me   /   tail -f render_all_<jobid>.log
# ============================================================================

#SBATCH --job-name=cmm523_render
#SBATCH --account=cusanovichlab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=10:00:00
#SBATCH --output=render_all_%j.log

SIF=/groups/darrenc/cmm_523/containers/cmm523_full.sif

cd "${SLURM_SUBMIT_DIR}"

TUTORIALS=(
  01_seurat_intro
  02_singler
  03_integration
  04_harmony
  05_pseudobulk
  06_trajectory
  07_signac
  08_rna_atac
)

echo "=========================================="
echo "Node:    $(hostname)"
echo "Started: $(date)"
echo "=========================================="

failed=()

for t in "${TUTORIALS[@]}"; do
    echo
    echo "########## ${t} ##########"
    echo "start: $(date +%H:%M:%S)"

    if apptainer exec -B /xdisk -B /groups "${SIF}" quarto render "${t}.qmd"; then
        echo "end:   $(date +%H:%M:%S)  OK"
    else
        echo "end:   $(date +%H:%M:%S)  FAILED"
        failed+=("${t}")
    fi
done

echo
echo "=========================================="
echo "Finished: $(date)"
if [ ${#failed[@]} -eq 0 ]; then
    echo "All tutorials rendered."
else
    echo "FAILED: ${failed[*]}"
fi
echo
echo "Rendered output:"
ls -la *.md
echo "=========================================="
