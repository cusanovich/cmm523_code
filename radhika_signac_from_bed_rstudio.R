dyn.load("/opt/ohpc/pub/apps/glpk/5.0/lib/libglpk.so.40")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")
dyn.load("/opt/ohpc/pub/apps/libpng/1.6.37/lib/libpng16.so.16")
dyn.load("/usr/lib64/atlas/libsatlas.so.3")

source("/opt/ohpc/admin/lmod/lmod/init/R")
Sys.setenv(MODULEPATH="/opt/ohpc/pub/moduledeps/gnu8-openmpi3:/opt/ohpc/pub/moduledeps/gnu8:/opt/ohpc/pub/modulefiles")
module("load python/3.11/3.11.4")

library(reticulate)
use_virtualenv('~/Python/pyenv_python_v3.11/')
library(Signac)
tester = CountFragments("/xdisk/darrenc/cmm_523/darrenc/radhika_proj/GSM5609895_31_fragments.tsv.gz")
cells = tester$CB[which(tester$frequency_count > 299)]
testobj = CreateFragmentObject(path = "/xdisk/darrenc/cmm_523/darrenc/radhika_proj/GSM5609895_31_fragments.tsv.gz",
  cells = cells)
peaks = CallPeaks(object=testobj,
                  macs2.path="~/Python/pyenv_python_v3.11/bin/macs3",
                  outdir ="/xdisk/darrenc/cmm_523/darrenc/radhika_proj")
updated = FeatureMatrix(fragments=testobj,features=peaks,cells=cells)
final_signac = CreateChromatinAssay(updated,fragments="/xdisk/darrenc/cmm_523/darrenc/radhika_proj/GSM5609895_31_fragments.tsv.gz")
saveRDS(final_signac,"/xdisk/darrenc/cmm_523/darrenc/radhika_proj/GSM5609895.rds")
