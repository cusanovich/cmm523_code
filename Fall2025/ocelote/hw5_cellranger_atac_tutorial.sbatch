#!/bin/bash
#SBATCH --job-name=atac_count
#SBATCH --account=darrenc
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00

cd /xdisk/darrenc/cmm_523/references/
wget https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz -O /xdisk/darrenc/cmm_523/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
tar -zxvf /xdisk/darrenc/cmm_523/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
rm /xdisk/darrenc/cmm_523/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz

wget https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_1k_nextgem/atac_pbmc_1k_nextgem_fastqs.tar  -O /xdisk/darrenc/cmm_523/references/atac_pbmc_1k_nextgem_fastqs.tar
tar -xvf /xdisk/darrenc/cmm_523/references/atac_pbmc_1k_nextgem_fastqs.tar
rm /xdisk/darrenc/cmm_523/references/atac_pbmc_1k_nextgem_fastqs.tar

mkdir /xdisk/darrenc/cmm_523/darrenc/cellranger_atac_count_tutorial/
cd /xdisk/darrenc/cmm_523/darrenc/cellranger_atac_count_tutorial/
cellranger-atac count --id=run_atac_count_pbmc_1k_nextgem \
                        --reference=/xdisk/darrenc/cmm_523/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                        --fastqs=/xdisk/darrenc/cmm_523/references/atac_pbmc_1k_nextgem_fastqs \
                        --sample=atac_pbmc_1k_nextgem