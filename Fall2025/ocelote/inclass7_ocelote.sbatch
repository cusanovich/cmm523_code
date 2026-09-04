#!/bin/bash
#SBATCH --job-name=inclass7
#SBATCH --account=darrenc
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00

module load R/4.2.2 gdal/3.3.2 geos/3.9.1 hdf5/1.10.5 libpng/1.6.37 glpk/5.0 libgit2/1.8.1

Rscript ~/ccsgab_projects/cmm523_code/Fall2025/ocelote/inclass7_ocelote.R