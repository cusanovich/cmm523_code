cd /xdisk/darrenc/cmm_523/darrenc/
mkdir run_cellranger_count
cd run_cellranger_count/

wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar

tar -xvf pbmc_1k_v3_fastqs.tar

cellranger count --help

cellranger count --id=run_count_1kpbmcs \
   --fastqs=/xdisk/darrenc/cmm_523/darrenc/run_cellranger_count/pbmc_1k_v3_fastqs \
   --sample=pbmc_1k_v3 \
   --transcriptome=/xdisk/darrenc/cmm_523/references/refdata-gex-GRCh38-2020-A

ls -1 run_count_1kpbmcs/outs
