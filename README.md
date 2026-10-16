# CMM 523 — Crash Course in Single-cell Genomics Analysis for Biologists

Course code, container recipes, and tutorials.

**Status:** rebuilt and verified October 2026. Not taught in 2026 — enrollment
was too low and the course was cancelled, so the material has been through a
full clean render but not through a classroom. The pinned versions date from
August 2026 and will need bumping before the next run. See
[Before the next run](#before-the-next-run).

---

## How the course is built

Four decisions shape everything else. They were made to stop the annual
rewrite, so they are worth keeping.

**Students build containers, not package installs.** Each in-class assignment
adds one Apptainer layer on top of the previous one, starting from an
instructor-built base image. Students still learn to install packages — from
CRAN, Bioconductor, and GitHub, including the token dance — but they do it
inside a recipe that is reproducible and that they can rebuild. The base image
holds the fiddly infrastructure that has nothing to teach.

**Versions are pinned, not floating.** The base image pins a Bioconductor
release and a Posit Package Manager snapshot date. The course does not drift
underneath itself mid-semester, and last year's tutorials still run.

**The tutorials are ours.** They are Quarto documents rendered to markdown and
read on GitHub. Rendering executes the code, so every figure in a `.md` is one
the code actually produced — if an upstream API breaks, the render fails and we
find out before the students do. Pointers to the original upstream tutorials are
kept in each one so students still learn to read other people's docs.

**Nothing caches.** Every tutorial sets `cache: false`. See
[Known fragilities](#known-fragilities) for why; this one was learned the hard
way and should not be undone to save a few minutes.

## Assignments, layers, and tutorials

Seven in-class assignments, seven homeworks, eight lectures, across ten sessions
in five weeks. Only five assignments add a container layer.

| Assignment | Image built | Adds | Tutorials |
|---|---|---|---|
| 1 | — (base image) | — | intro, HPC, containers |
| 2 | `assignment2.sif` | Seurat, SeuratObject | 01 Seurat intro |
| 3 | `assignment3.sif` | harmony, SingleR, celldex, scRNAseq, scrapper, SeuratData | 02 SingleR, 03 integration, 04 harmony |
| 4 | `assignment4.sif` | edgeR, statmod | 05 pseudobulk, 06 trajectory |
| 5 | — (reuses `assignment4.sif`) | — | 06 trajectory |
| 6 | `assignment6.sif` | Signac, EnsDb.Hsapiens.v86, biovizBase, hdf5r | 07 Signac, 08 RNA+ATAC |
| 7 | `assignment7.sif` | chromVAR, motifmatchr, TFBSTools, JASPAR2024, BSgenome hg38, presto | 09 WNN |

Each layer bootstraps from the one before it, so `assignment7.sif` accumulates
everything. That is the image students keep.

Assignment 3 is where GitHub tokens are taught — it is the first layer that
installs from GitHub, and the shared cluster IP hits GitHub's unauthenticated
rate limit (60/hr) almost immediately with a class of twenty.

monocle3 is pre-baked into the base image rather than being an assignment layer.
It does not import Seurat, so baking it in takes nothing away from assignment 2,
and building it from source is slow and failure-prone.

## Repository layout

```
container/
  cmm523_base.def                 base image: R, Bioconductor, Quarto, monocle3
  cmm523_full.def                 base + everything students install
  student/assignment*.def         the five student layers
  build_layer.sbatch              build any recipe as a batch job
  fetch_references.sbatch         mirror course datasets to /groups
  rserver_cmm523.sh               RStudio Server from a container, in OOD Desktop
tutorials/
  NN_name.qmd                     source
  NN_name.md                      rendered — this is what students read
  figs/                           generated figures
  seurat_to_cds.R                 Seurat -> monocle3 helper (replaces SeuratWrappers)
  render_tutorial.sbatch          render one
  render_all.sbatch               render all nine, in order
Fall2025/                         archived prior version
```

## Storage on the cluster

```
/groups/darrenc/cmm_523/
  containers/    cmm523_base.sif, cmm523_full.sif
  recipes/       staged copies students read (keep world-readable)
  references/    mirrored datasets — panc8, pbmcMultiome, hcabm40k,
                 10x multiome, Cell Ranger references, edgeR SeuratObj.rds
/xdisk/darrenc/cmm_523/<netid>/    student build and work area
```

`/groups` for anything durable. `/xdisk` is wiped roughly every ten months, so
nothing that takes an hour to rebuild should live there.

Access depends on students being in the `cusanovichlab` **Unix group**. Nothing
works without it.

## Running things

```bash
# interactive R in a container
apptainer exec -B /xdisk -B /groups $SIF R

# RStudio, from an OOD Desktop session (not the OOD RStudio app —
# it cannot launch containers)
bash rserver_cmm523.sh /path/to/your.sif

# render one tutorial
sbatch render_tutorial.sbatch tutorials/07_signac.qmd

# render all nine
sbatch render_all.sbatch
```

Memory on Puma is allocated **per core** (4 GB each) and there is no `--mem`
flag. To get 32 GB, ask for 8 cores.

### Render times

Measured on a clean full run, 30 September 2026, 16 cores. Nothing is cached,
so these are the real cost every time.

| Tutorial | Time | | Tutorial | Time |
|---|---|---|---|---|
| 01 seurat_intro | 1 min | | 06 trajectory | 7 min |
| 02 singler | 4 min | | 07 signac | 29 min |
| 03 integration | 5 min | | 08 rna_atac | 19 min |
| 04 harmony | 5 min | | 09 wnn | 90 min |
| 05 pseudobulk | 2 min | | **all nine** | **~2h 45m** |

`render_all.sbatch` requests 16 cores and 10 hours, which is ample. Note that
these timings vary by 30–40% between nodes, so treat them as a guide rather
than a budget.

## Before the next run

Steps 1–5 are the annual version bump; the rest is verification.

1. **Pick the Bioconductor release.** Releases land in April and October. Pin
   whichever is current when the semester starts so it cannot change
   mid-course. Confirm the matching `bioconductor/bioconductor_docker:RELEASE_X_Y`
   tag exists.

2. **Update `cmm523_base.def`:**
   - the `From:` base image tag
   - the Posit snapshot date — pick one *after* the Bioconductor release
   - the Ubuntu codename in the PPM URL (`__linux__/noble/`) — it must match
     the base image's Ubuntu version, or you silently get source packages and
     the build goes from ~25 minutes to ~90
   - `QUARTO_VERSION`

3. **Build the base image.** `sbatch build_layer.sbatch cmm523_base.def cmm523_base.sif`,
   about 25 minutes. The recipe has a guard that fails the build if Seurat
   appears as a side effect — if it trips, something new started importing
   Seurat and assignment 2 needs rethinking.

4. **Build `cmm523_full.sif`.** Needed for rendering, and the fallback if a
   student's chain breaks mid-semester.

5. **Re-render all nine tutorials and read every one.** Budget three hours.
   This is the step that catches upstream API changes, and there is no
   substitute for reading the output — a deprecation warning renders as
   ordinary text and will not fail the job.

6. **Check for orphaned figures** before committing the re-render. Chunks that
   get renamed or removed leave their old `figs/` files behind, and they
   accumulate across years:

   ```bash
   cd tutorials
   grep -ho 'figs/[A-Za-z0-9._-]*' *.md | sort -u > /tmp/ref.txt
   ls figs/ | sed 's|^|figs/|' | sort > /tmp/have.txt
   comm -13 /tmp/ref.txt /tmp/have.txt    # exist but unreferenced -- delete
   comm -23 /tmp/ref.txt /tmp/have.txt    # referenced but missing -- a real problem
   ```

7. **Build the five student layers in order**, base → 2 → 3 → 4 → 6 → 7. This
   verifies the chain students actually walk, which the full image does not.

8. **Verify the reference mirror** is intact and re-run `fetch_references.sbatch`
   for anything missing.

9. **Re-stage recipes** to `/groups/darrenc/cmm_523/recipes/` and confirm they
   are world-readable, including every parent directory (`namei -l` is the
   quick check).

10. **Update the assignment documents** for any command that changed.

11. **Update the lecture decks** — dates on all eight title slides, and the
    journal club section number on Lecture 8 slide 29.

12. **Add students to `cusanovichlab`** once the roster is final.

## Known fragilities

Things that have broken before and will break again.

**The reference you choose decides the answer.** This is the biggest one, and
it is invisible. Tutorial 2 originally used `HumanPrimaryCellAtlasData`, which
contains a `DC` label — and assigned it to **zero** cells, because HPCA is bulk
microarray of sorted and cultured populations whose dendritic cell profile does
not resemble blood cDCs or pDCs. The dendritic cells were real and present;
they were simply called monocytes, with no warning, no pruning, and nothing in
the output to suggest a problem. Switching to `MonacoImmuneData`, which is
RNA-seq of sorted blood populations, labelled 2% of cells as dendritic cells —
the expected frequency. **Match the reference to the tissue**, and always check
the annotation against a population you know should be there.

**Caching restores objects without re-running the chunks that made them.**
knitr re-runs a chunk only when *that chunk's own code* changes; it has no idea
that an object the chunk consumes was rebuilt upstream. Editing the reference in
tutorial 2 therefore left a downstream chunk restoring a stale object, and the
rendered figure showed August's annotation beside October's in the same panel.
Nothing errored. Every tutorial now sets `cache: false`. Tutorials 7, 8 and 9
additionally *cannot* cache — knitr's lazy-load database will not hold an object
over 2 GB, and a ChromatinAssay over 100,000 peaks exceeds it.

**Packages that install clean and fail at runtime.** `R.utils`, `scrapper`,
`statmod`, `hdf5r`, `biovizBase`, `glmGamPoi`. All are Suggests rather than
Imports of something the course uses, so nothing complains until the call
happens. They are installed explicitly for this reason — do not "clean them up."

**`options(repos = ...)` replaces rather than amends.** Setting CRAN with
`c(CRAN = ...)` wipes the Bioconductor repos, after which Bioconductor packages
are silently skipped. Always amend: `local({r <- getOption("repos"); r["CRAN"] <- ...})`.

**`remotes::install_github()` resolves against CRAN only** unless you first set
`options(repos = BiocManager::repositories())`. Any GitHub package with
Bioconductor imports fails without it.

**Genome build mismatches fail silently.** `EnsDb.Hsapiens.v86` is hg38, `v75`
is hg19. Pairing hg19 data with hg38 annotation produces TSS enrichment near
zero, the QC filter then removes every cell, and the error surfaces much later
as something unrelated about UMAP neighbours. This is now a teaching point in
tutorial 7 and HW6.

**`~/.Renviron` shadows the container library.** A stale `R_LIBS` in a user's
home directory overrides the container's library path and produces
`undefined symbol` errors that look like a broken build. Neither `--env` nor
`Rprofile.site` reliably wins. Move the file aside.

**`de.method` has to match the reference type.** SingleR's default `"classic"`
is built for bulk references, which both HPCA and Monaco are. `"wilcox"` and
`"t"` are for single-cell references with hundreds of cells per label. Asking
for a Wilcoxon test on four sorted samples gives unstable marker selection, and
since markers drive the correlation scores, that becomes confident wrong labels.

**API changes to watch:** Seurat 5 uses `layer=`, not `slot=` (defunct, errors
rather than warns). `RunChromVAR` was removed in Signac 1.17 — call chromVAR
directly. JASPAR bumps versions and `getMatrixSet()` needs an RSQLite
connection, not the package object.

**chromVAR under BiocParallel runs out of memory** and reports it as
`wrong args for environment subassignment`. Register `SerialParam()`.

## Open items

- **`cellranger-atac` has no module on Puma.** `cellranger/9.0.1` exists;
  the ATAC counterpart does not. InClass5 is written to work either way with an
  instructor note. Resolve by requesting a module from HPC or staging a shared
  install under `/groups/darrenc/cmm_523/tools/`.
- **`/groups/darrenc` is group-writable with no sticky bit.** Once students are
  in `cusanovichlab` they can rename or delete top-level entries there.
  `chmod +t /groups/darrenc` fixes it. Also worth checking what else in the lab
  directory becomes readable to the class.
- **The 2025 scripts are still in the repository root** — `WNN_tutorial.R`,
  `Signac_tutorial.R`, `edgeR_tutorial.R`, and the rest. They contain
  `RunChromVAR`, `JASPAR2020`, and `slot=`, all of which now error. A student
  who finds them first will get failures that have nothing to do with the
  course. Move them into `Fall2025/` or delete them.
- **A pre-flight check, like the one CMM 522 has.** Less urgent here because the
  container pins most of what would otherwise drift, but a script that confirms
  the base image still builds, the layer chain still bootstraps, and the
  reference mirror is intact would turn steps 3–8 above into one job.
