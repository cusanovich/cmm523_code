#!/bin/bash
# ============================================================================
# CMM 523 -- pre-flight check
#
# Run this in August, before touching anything:
#
#     cd ~/ccsgab_projects/cmm523_code/container
#     sbatch preflight_cmm523.sbatch
#
# It answers one question: is last year's setup still intact? Run it BEFORE
# the version bump in the README, so that when something breaks during the
# bump you know it was the bump that broke it.
#
# It is read-only. It builds nothing, downloads nothing, and changes no file.
# Expect two to five minutes, nearly all of it waiting on the link checks.
#
# What it does NOT do, because each takes hours: rebuild the base image,
# build the five student layers, or re-render the tutorials. Those are the
# annual work. This confirms the ground you are standing on first.
# ============================================================================

set -u
FAIL=0
note()  { printf '\n\033[1m== %s ==\033[0m\n' "$1"; }
pass()  { printf '  ok    %s\n' "$1"; }
warn()  { printf '  WARN  %s\n' "$1"; FAIL=1; }
bad()   { printf '  FAIL  %s\n' "$1"; FAIL=1; }

# --- what the recipes currently pin -------------------------------------
# Update these when you deliberately move the course forward.
BIOC_TAG="RELEASE_3_23"
PPM_SNAPSHOT="2026-08-01"
QUARTO_VERSION="1.7.32"

# NB: do NOT name this variable GROUPS. In bash that is a read-only special
# array of the current user's group IDs; assigning to it is silently ignored
# and ${GROUPS} then expands to your primary GID, so every path built from it
# becomes garbage. Learned the hard way.
COURSE=/groups/darrenc/cmm_523
CONTAINERS="${COURSE}/containers"
RECIPES="${COURSE}/recipes"
REFS="${COURSE}/references"

# Print the resolved paths. If a variable ever fails to expand the way it
# looks like it should, this line makes it obvious instead of producing a
# page of confusing FAILs.
note "Paths"
printf '  %s\n' "$CONTAINERS" "$RECIPES" "$REFS"
if [ ! -d "$COURSE" ]; then
  bad "course directory does not exist: $COURSE"
  echo "  Everything below would fail for that reason alone. Stopping."
  exit 1
fi
pass "course directory present"

note "Apptainer"
if command -v apptainer >/dev/null 2>&1; then
  pass "$(apptainer --version)"
else
  bad "apptainer not found -- are you on a compute node?"
  exit 1
fi

# --- container images -----------------------------------------------------
note "Container images"
for img in cmm523_base.sif cmm523_full.sif; do
  p="${CONTAINERS}/${img}"
  if [ -r "$p" ]; then
    pass "$(ls -lh "$p" | awk '{print $5, $9}')"
  else
    bad "missing or unreadable: $p"
  fi
done

# --- the images actually work --------------------------------------------
note "Base image contents"
if [ -r "${CONTAINERS}/cmm523_base.sif" ]; then
  apptainer exec -B /groups "${CONTAINERS}/cmm523_base.sif" Rscript -e '
    cat("        ", R.version.string, "\n", sep = "")
    cat("        Bioconductor ", as.character(BiocManager::version()), "\n", sep = "")
    cat("        packages installed: ", nrow(installed.packages()), "\n", sep = "")
    for (p in c("monocle3", "BiocManager", "knitr", "rmarkdown")) {
      cat(if (requireNamespace(p, quietly = TRUE)) "  ok    " else "  FAIL  ", p, "\n", sep = "")
    }
    # The base image must NOT contain Seurat -- installing it is assignment 2.
    if (requireNamespace("Seurat", quietly = TRUE)) {
      cat("  FAIL  Seurat is present in the BASE image; assignment 2 has nothing to do\n")
    } else {
      cat("  ok    Seurat absent from base, as intended\n")
    }
  ' 2>&1 | grep -E '^( |  )' || bad "could not run R inside the base image"
  apptainer exec -B /groups "${CONTAINERS}/cmm523_base.sif" \
    Rscript -e 'quit(status = !requireNamespace("monocle3", quietly=TRUE))' || FAIL=1
fi

note "Full image contents"
if [ -r "${CONTAINERS}/cmm523_full.sif" ]; then
  apptainer exec -B /groups "${CONTAINERS}/cmm523_full.sif" Rscript -e '
    pkgs <- c("Seurat","SeuratObject","harmony","SingleR","celldex","scRNAseq",
              "scrapper","edgeR","statmod","Signac","EnsDb.Hsapiens.v86",
              "biovizBase","hdf5r","chromVAR","motifmatchr","TFBSTools",
              "JASPAR2024","BSgenome.Hsapiens.UCSC.hg38","SeuratData","presto",
              "monocle3","knitr","rmarkdown")
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) {
      cat("  FAIL  missing: ", paste(missing, collapse = ", "), "\n", sep = "")
      quit(status = 1)
    }
    cat("  ok    all ", length(pkgs), " course packages load\n", sep = "")
  ' || FAIL=1
fi

# --- recipes staged and readable -----------------------------------------
note "Staged recipes"
EXPECTED_RECIPES="assignment2_seurat.def assignment3_annotation.def
assignment4_pseudobulk.def assignment6_signac.def assignment7_motifs.def
cmm523_base.def cmm523_full.def rserver_cmm523.sh"
for r in $EXPECTED_RECIPES; do
  if [ -r "${RECIPES}/${r}" ]; then pass "$r"; else bad "missing from ${RECIPES}: $r"; fi
done

note "Readable by students"
echo "  every component needs x for others, and the file needs r"
if [ -e "${RECIPES}/cmm523_base.def" ]; then
  namei -l "${RECIPES}/cmm523_base.def" | sed 's/^/        /'
fi
# /groups/darrenc is group-writable; without the sticky bit anyone in
# cusanovichlab can rename or delete top-level entries there.
PERM=$(stat -c '%A' /groups/darrenc 2>/dev/null)
echo "        /groups/darrenc  $PERM"
case "$PERM" in
  *t|*T) pass "sticky bit set -- students cannot rename or remove others' entries" ;;
  *)     warn "no sticky bit on /groups/darrenc -- run: chmod +t /groups/darrenc" ;;
esac

# --- reference mirror -----------------------------------------------------
note "Reference mirror"
# Existence is not enough. A failed or interrupted download leaves an empty
# directory behind, which looks fine to -e and is useless to a tutorial. Check
# that something is actually in there.
for d in refdata-gex-GRCh38-2024-A refdata-cellranger-arc-GRCh38-2024-A \
         Rdatalib edgeR_pal pbmc_multiome; do
  t="${REFS}/${d}"
  if [ ! -e "$t" ]; then
    bad "missing: $t   (re-run fetch_references.sbatch)"
    continue
  fi
  # Some of these are symlinks into the lab's shared reference collection
  # rather than copies. du and find must be told to follow them (-L), or a
  # perfectly good symlink reports as a few KB containing no files.
  link=""
  if [ -L "$t" ]; then
    tgt=$(readlink -f "$t" 2>/dev/null)
    if [ -e "$tgt" ]; then
      link="  -> $tgt"
    else
      bad "BROKEN LINK: $t -> $(readlink "$t")"
      continue
    fi
  fi
  size=$(du -shL "$t" 2>/dev/null | awk '{print $1}')
  kb=$(du -skL  "$t" 2>/dev/null | awk '{print $1}')
  nfiles=$(find -L "$t" -type f 2>/dev/null | head -2 | wc -l)
  if [ "$nfiles" -eq 0 ]; then
    bad "EMPTY: $t holds no files ($size)${link}"
    echo "        the download never completed -- re-run fetch_references.sbatch"
  elif [ "${kb:-0}" -lt 1024 ]; then
    warn "$size  $d   -- suspiciously small${link}"
  else
    pass "$size  $d${link}"
  fi
done

# --- upstream pins still exist -------------------------------------------
note "Upstream pins"
code=$(curl -s -o /dev/null -w '%{http_code}' --max-time 20 \
  "https://hub.docker.com/v2/repositories/bioconductor/bioconductor_docker/tags/${BIOC_TAG}")
case "$code" in
  200) pass "bioconductor_docker:${BIOC_TAG} still published" ;;
  *)   warn "bioconductor_docker:${BIOC_TAG} returned $code -- confirm the tag by hand" ;;
esac

code=$(curl -s -o /dev/null -w '%{http_code}' -L --max-time 20 \
  "https://packagemanager.posit.co/cran/__linux__/noble/${PPM_SNAPSHOT}/src/contrib/PACKAGES.gz")
case "$code" in
  2*|30*) pass "Posit snapshot ${PPM_SNAPSHOT} still served" ;;
  *)      warn "Posit snapshot ${PPM_SNAPSHOT} returned $code" ;;
esac

code=$(curl -s -o /dev/null -w '%{http_code}' -L --max-time 20 \
  "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb")
case "$code" in
  2*|30*) pass "Quarto ${QUARTO_VERSION} release still downloadable" ;;
  *)      warn "Quarto ${QUARTO_VERSION} returned $code" ;;
esac

# --- data the tutorials fetch at render time ------------------------------
note "Data and links the tutorials use"
echo "  (a redirect is fine; 404 and 5xx are not)"
while IFS= read -r u; do
  [ -z "$u" ] && continue
  code=$(curl -s -o /dev/null -w '%{http_code}' -L --max-time 25 "$u" 2>/dev/null)
  # Several of these URLs share a host, and hammering one host in quick
  # succession gets the checker throttled -- which looks exactly like a dead
  # link. On no response, wait and try once more before believing it.
  if [ "$code" = "000" ]; then
    sleep 5
    code=$(curl -s -o /dev/null -w '%{http_code}' -L --max-time 25 "$u" 2>/dev/null)
  fi
  case "$code" in
    2*|30*)      pass "$code  $u" ;;
    401|403|405) pass "$code  $u   (blocks automated checks; open it by hand)" ;;
    000)         warn "no response after retry  $u"
                 printf '        a compute node cannot always reach everything; try a browser\n' ;;
    *)           bad  "$code  $u" ;;
  esac
  sleep 1
done <<'URLS'
https://bioinf.wehi.edu.au/edgeR/UserGuideData/SeuratObj.rds
https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2024-A.tar.gz
https://bioconductor.org/books/release/SingleRBook/
https://bioconductor.org/packages/release/bioc/html/edgeR.html
https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
https://portals.broadinstitute.org/harmony/
https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
https://satijalab.org/seurat/articles/integration_introduction.html
https://satijalab.org/seurat/articles/integration_mapping.html
https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis
https://stuartlab.org/signac/articles/pbmc_vignette
https://stuartlab.org/signac/articles/pbmc_multiomic
https://stuartlab.org/signac/articles/motif_vignette
URLS

# --- summary --------------------------------------------------------------
note "Summary"
if [ "$FAIL" -eq 0 ]; then
  echo "  Last year's setup is intact. Safe to start the version bump."
else
  echo "  Something above needs attention. Fix it BEFORE bumping versions,"
  echo "  so you are not debugging two problems at once."
fi
exit "$FAIL"
