#!/bin/bash
#SBATCH --job-name=mkfastq
#SBATCH --account=darrenc
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00

module load bcl2fastq2

#cd ~/bin
#curl -o cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1758972007&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=a6uFr5IreOsY0DXGmvFuSzljljF8gMW9CiMjX08XcXtn4ne21Ll9WaZb50~cA0YwxTtYohGrW5qamoXfPF7~wUJ0W0cn7XHdeWWjpQjP9GPUqp0WJ0XdpqCov1Tqucd1qaeBcVtau41285DfaIvPwh1CTmRqzb-lFk4WDUuCBPUF7L950-ayAHtkdGW9qlG4FfmG76YiVtywUgo9Z~HICW0XzZzC-QYdsYZylyyu0QGM8RGAS4p6cv4sKQ8g5DB36MFlOaNZV6fkFM5aX4iKHJzGLh8X9SgnYtOL~5ls3rUW~0~ljDj1s~GH7Ynkfh3YXjC1yXcxofmG-CGH~qyemg__"

#tar -xzvf cellranger-9.0.1.tar.gz

cd /xdisk/darrenc/cmm_523/darrenc/
mkdir cellranger_mkfastq_tutorial
cd cellranger_mkfastq_tutorial/

wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv
tar -zxvf cellranger-tiny-bcl-1.2.0.tar.gz
rm cellranger-tiny-bcl-1.2.0.tar.gz

cat cellranger-tiny-bcl-simple-1.2.0.csv
tree -L 2 cellranger-tiny-bcl-1.2.0/
cellranger mkfastq --help

 cellranger mkfastq --id=tutorial_walk_through \
   --run=/xdisk/darrenc/cmm_523/darrenc/cellranger_mkfastq_tutorial/cellranger-tiny-bcl-1.2.0 \
   --csv=/xdisk/darrenc/cmm_523/darrenc/cellranger_mkfastq_tutorial/cellranger-tiny-bcl-simple-1.2.0.csv

cd /xdisk/darrenc/cmm_523/darrenc/cellranger_mkfastq_tutorial/tutorial_walk_through/outs/fastq_path
ls -1

echo "Done!"