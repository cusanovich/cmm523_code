#!/bin/bash
# ============================================================================
# CMM 523 -- run the pre-flight check as a batch job
#
#     cd ~/ccsgab_projects/cmm523_code/container
#     sbatch preflight_cmm523.sbatch
#
# Result lands in preflight_<jobid>.log
# ============================================================================

#SBATCH --job-name=cmm523_preflight
#SBATCH --account=cusanovichlab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=00:30:00
#SBATCH --output=preflight_%j.log

cd "${SLURM_SUBMIT_DIR}"

echo "Submitted from: ${SLURM_SUBMIT_DIR}"
echo "Running on:     $(hostname)"
echo "Started:        $(date)"
echo

bash preflight_cmm523.sh
STATUS=$?

echo
echo "Finished: $(date)"
if [ "$STATUS" -eq 0 ]; then
    echo "Pre-flight passed. Last year's setup is intact."
else
    echo "Pre-flight found problems. Search this log for 'WARN' and 'FAIL'."
fi
exit "$STATUS"
