module load bcl2fastq2

cd /xdisk/darrenc/cmm_523/darrenc/
mkdir run_cellranger_mkfastq
cd run_cellranger_mkfastq

wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv
tar -zxvf cellranger-tiny-bcl-1.2.0.tar.gz

cat cellranger-tiny-bcl-simple-1.2.0.csv
tree -L 2 cellranger-tiny-bcl-1.2.0/

cellranger mkfastq --help

cellranger mkfastq --id=tutorial_walk_through \
  --run=/xdisk/darrenc/cmm_523/darrenc/run_cellranger_mkfastq/cellranger-tiny-bcl-1.2.0 \
  --csv=/xdisk/darrenc/cmm_523/darrenc/run_cellranger_mkfastq/cellranger-tiny-bcl-simple-1.2.0.csv

cd /xdisk/darrenc/cmm_523/darrenc/run_cellranger_mkfastq/tutorial_walk_through/outs/fastq_path
ls -1

ls -1 H35KCBCXY/test_sample