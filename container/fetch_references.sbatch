#!/bin/bash
# ============================================================================
# CMM 523 -- stage the shared reference genomes
#
#   sbatch fetch_references.sbatch
#
# Downloads and unpacks the Cell Ranger RNA and ATAC/ARC references into the
# shared course directory, then makes them world-readable so students can use
# them without their own copies.
#
# WHY THIS EXISTS: /xdisk allocations are purged periodically, so these
# disappear between semesters. Run this once before the course starts. Roughly
# 21 GB of download and about 40 GB unpacked -- check `uquota` first.
#
# These URLs are 10x's public CDN and do not expire, unlike the signed links
# from the download page that require the End User Agreement form.
# ============================================================================

#SBATCH --job-name=cmm523_refs
#SBATCH --account=cusanovichlab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=fetch_refs_%j.log

set -e

REFDIR=/groups/darrenc/cmm_523/references
mkdir -p "${REFDIR}"
cd "${REFDIR}"

echo "Started: $(date)"
echo "Target:  ${REFDIR}"
df -h "${REFDIR}" | tail -1

# --- RNA reference (used by cellranger count) -------------------------------
# Check https://www.10xgenomics.com/support/software/cell-ranger/downloads for
# the current version before running; bump this deliberately once a year.
RNA_REF=refdata-gex-GRCh38-2024-A

if [ -d "${RNA_REF}" ]; then
    echo "Already have ${RNA_REF}, skipping."
else
    echo "=== downloading ${RNA_REF} (~11 GB) ==="
    wget -q --show-progress "https://cf.10xgenomics.com/supp/cell-exp/${RNA_REF}.tar.gz"
    tar -xzf "${RNA_REF}.tar.gz"
    rm "${RNA_REF}.tar.gz"
fi

# --- ATAC / multiome reference (used by cellranger-atac count) --------------
# This is the ARC reference, which serves both cellranger-atac and
# cellranger-arc. It is a different file from the RNA reference above --
# the two pipelines cannot share one. Note the URL path is cell-arc, not
# cell-atac, and the naming convention differs from the GEX reference.
ATAC_REF=refdata-cellranger-arc-GRCh38-2024-A

if [ -d "${ATAC_REF}" ]; then
    echo "Already have ${ATAC_REF}, skipping."
else
    echo "=== downloading ${ATAC_REF} (~10 GB) ==="
    wget -q --show-progress "https://cf.10xgenomics.com/supp/cell-arc/${ATAC_REF}.tar.gz"
    tar -xzf "${ATAC_REF}.tar.gz"
    rm "${ATAC_REF}.tar.gz"
fi

# --- make readable by the class --------------------------------------------
chmod -R a+rX "${REFDIR}"

echo
echo "=== done ==="
du -sh "${REFDIR}"/*
echo "Finished: $(date)"
