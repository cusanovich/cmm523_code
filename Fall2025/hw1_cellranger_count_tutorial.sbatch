#!/bin/bash
#SBATCH --job-name=count
#SBATCH --account=darrenc
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00

#wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar -O /xdisk/darrenc/cmm_523/references/pbmc_1k_v3_fastqs.tar
#tar -xvf /xdisk/darrenc/cmm_523/references/pbmc_1k_v3_fastqs.tar
#rm /xdisk/darrenc/cmm_523/references/pbmc_1k_v3_fastqs.tar

#wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz -O /xdisk/darrenc/cmm_523/references/refdata-gex-GRCh38-2020-A.tar.gz
#tar -zxvf /xdisk/darrenc/cmm_523/references/refdata-gex-GRCh38-2020-A.tar.gz
#rm /xdisk/darrenc/cmm_523/references/refdata-gex-GRCh38-2020-A.tar.gz

cd /xdisk/darrenc/cmm_523/darrenc/cellranger_count_tutorial/
cellranger count --id=run_count_1kpbmcs \
   --fastqs=/xdisk/darrenc/cmm_523/references/pbmc_1k_v3_fastqs \
   --sample=pbmc_1k_v3 \
   --transcriptome=/xdisk/darrenc/cmm_523/references/refdata-gex-GRCh38-2020-A \
   --create-bam true

ls -1 run_count_1kpbmcs/outs