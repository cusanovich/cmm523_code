#!/bin/bash
#SBATCH --job-name=signac
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=darrenc

module load python/3.11 R/4.2.2 gdal/3.3.2 geos/3.9.1 hdf5/1.10.5 libpng/1.6.37 glpk/5.0 libgit2/1.8.1
source ~/Python/pyenv_python_v3.11/bin/activate
Rscript ~/ccsgab_projects/cmm523_code/radhika_signac_from_bed.R