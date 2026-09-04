#!/bin/bash
#SBATCH --job-name=atac_count
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=4gb
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=darrenc

cd /xdisk/darrenc/cmm_523/darrenc/cellranger_atac_tutorial/
#wget https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_1k_nextgem/atac_pbmc_1k_nextgem_fastqs.tar -O /xdisk/darrenc/cmm_523/darrenc/cellranger_atac_tutorial/atac_pbmc_1k_nextgem_fastqs.tar
#tar -xvf /xdisk/darrenc/cmm_523/darrenc/cellranger_atac_tutorial/atac_pbmc_1k_nextgem_fastqs.tar
cellranger-atac count --id=atac_pbmc_1k_nextgem \
                        --reference=/xdisk/darrenc/cmm_523/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                        --fastqs=/xdisk/darrenc/cmm_523/darrenc/cellranger_atac_tutorial/atac_pbmc_1k_nextgem_fastqs/ \
                        --sample=atac_pbmc_1k_nextgem