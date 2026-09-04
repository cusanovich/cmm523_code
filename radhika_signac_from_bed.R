library(Signac)
tester = CountFragments("/xdisk/darrenc/cmm_523/darrenc/radhika_proj/GSM5609895_31_fragments.tsv.gz")
cells = tester$CB[which(tester$frequency_count > 299)]
testobj = CreateFragmentObject(path = "/xdisk/darrenc/cmm_523/darrenc/radhika_proj/GSM5609895_31_fragments.tsv.gz",
  cells = cells)
peaks = CallPeaks(object=testobj,
                  macs2.path="~/Python/pyenv_python_v3.11/bin/macs3",
                  outdir ="/xdisk/darrenc/cmm_523/darrenc/radhika_proj/")
updated = FeatureMatrix(fragments=testobj,features=peaks,cells=cells)
final_signac = CreateChromatinAssay(updated,fragments="/xdisk/darrenc/cmm_523/darrenc/radhika_proj/GSM5609895_31_fragments.tsv.gz")
saveRDS(final_signac,"/xdisk/darrenc/cmm_523/darrenc/radhika_proj/GSM5609895.rds")
