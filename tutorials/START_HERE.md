# CMM 523 — working through it on your own

This course is normally taught live, over ten sessions in five weeks, with an
instructor in the room. You are doing it self-paced, which mostly works, but a
few things are written assuming someone is there to ask. This page covers those.

---

Your HPC account and `cusanovichlab` membership are already set up. If
`ls /groups/darrenc/cmm_523/recipes/` ever says permission denied, something
changed — tell Darren.

## Where everything lives

| What | Where |
|---|---|
| Tutorials (read these) | `github.com/cusanovich/cmm523_code` → `tutorials/*.md` |
| Container recipes | `/groups/darrenc/cmm_523/recipes/` |
| Instructor images | `/groups/darrenc/cmm_523/containers/` |
| Reference data | `/groups/darrenc/cmm_523/references/` |
| Your working area | `/xdisk/darrenc/cmm_523/<your_netid>/` |

The tutorials render on GitHub with all their figures, so you can read them in a
browser without running anything. **Read the tutorial before the assignment that
uses it.** The assignments assume you have.

## Day one: build the whole container chain

Do this before any analysis.

Each in-class assignment has you add one layer to an Apptainer container,
building on the previous one. The chain is:

```
cmm523_base.sif  →  assignment2  →  assignment3  →  assignment4  →  assignment6  →  assignment7
   (provided)        Seurat         annotation      edgeR           Signac          motifs
```

Because each layer depends on the one before it, a failure at assignment 2
blocks everything after it. In a live class you would find that out with the
instructor standing there. On your own, you want to find out on day one.

So: work through the container-building steps of assignments 2, 3, 4, 6 and 7
in order, *before* doing any of the analysis. Submit each as a batch job with
`build_layer.sbatch`. If all five images build, the rest of the course is
unblocked and you can go back to assignment 1 and work through properly.

Two things to know before you start:

- **Assignment 3 needs a GitHub personal access token.** The assignment explains
  how to make one. Without it the build fails partway through, which looks like
  a container problem and is not.
- **Budget around 15 GB on `/xdisk`.** Each layer is a complete image, not a
  diff, so you end up with five images of roughly 3 GB each. Check your quota
  before you start, and tell Darren what the real total comes to — that number
  matters for running this with a full class.

## Then work through in order

Seven in-class assignments and seven homeworks. Do them in order: in-class *N*,
then homework *N*, then in-class *N+1*. The homeworks extend the in-class work
and sometimes depend on files it produced.

Assignment 5 adds no container layer — it reuses `assignment4.sif`.

The eight lecture decks give the conceptual background for each session. They
were written to be presented rather than read, so expect them to be sparse.
Skim the deck for a session before the assignment, mainly to see what the
session is about and what vocabulary to expect.

## When you get stuck

The assignments sometimes say "ask for help" or "we will check in here." You do
not have that, so:

**Read the error from the bottom up.** The last line is usually the symptom; the
cause is often several lines above it.

**Check the tutorial.** Most of what an in-class assignment asks you to do is
demonstrated in the matching tutorial with working code.

**Suspect memory.** Memory on Puma comes with cores — 4 GB each, and there is no
`--mem` flag. If R dies with no error message at all and the session simply
stops, you ran out. Ask for more cores.

**Suspect a stale `~/.Renviron`.** If R inside a container reports missing
packages that you know are installed, or something about an `undefined symbol`,
a leftover `R_LIBS` in your home directory is probably overriding the
container's library. Move the file aside and try again.

**Then ask Darren** — but include what you tried and the actual error text, not
just "it didn't work."

## What to write down

This is the first time anyone has worked through this material without the
instructor present, so the most useful thing you can produce is a record of
where it was unclear.

Keep a running note. For each point where you stopped:

- which assignment and step
- what you expected to happen
- what happened instead
- what you tried, and what eventually worked

"I had to guess what this meant" is a useful entry. So is "this took an hour and
I expected ten minutes." Do not tidy it up or decide something is too small to
mention — the small confusions are the ones that are invisible from the
instructor's side, because he already knows the answer.
