# Prediction accuracy against RefSeq

This note records how Orphos's predictions compare with RefSeq CDS
annotations, where the errors are, and what has been tried to reduce them.
Orphos's output is byte-identical to Prodigal (GitHub master), so every number
below also describes Prodigal.

## Benchmark

- `scripts/refseq_cds_tsv.py` converts an NCBI feature table
  (`efetch ... rettype=ft`) into `tests/data/refseq/<accession>.cds.tsv`,
  keeping only complete, single-interval, non-pseudo CDS.
- `scripts/evaluate_accuracy.py` scores GFF predictions against those tables.
- `cargo test -p orphos-cli --test reference_annotations` runs the same
  comparison for every complete-genome fixture and snapshots the summary.

A reference gene is *found* when a prediction shares its 3' end (the stop
codon) and *exact* when the 5' end matches too. Predictions with no reference
gene at their 3' end are *false positives* relative to the annotation.

RefSeq annotations are themselves produced by gene finders plus protein
evidence and curation (EcoCyc for *E. coli*, SubtiWiki for *B. subtilis*),
so they are a strong proxy rather than ground truth: a disagreement is not
always an Orphos error.

## Results (34 genomes, single mode)

| | value |
| --- | ---: |
| reference CDS | 85,136 |
| found by 3' end | 97.6% |
| exact (start and stop) | 86.9% |
| start correct, given found | 89.0% |
| false positives / predictions | 7.5% |

Per-genome numbers for the fixtures shipped in the repository are in
`tests/snapshots/reference_annotations__refseq_cds_accuracy_summary.snap`.
Start accuracy ranges from 74% (*Treponema pallidum*) and 77%
(*Mycobacterium tuberculosis*) to 95% (*Haemophilus influenzae*); the
GC-rich genomes are the weakest.

## Where the errors are

- **Start choice dominates.** Of the reference genes that are found, 11% have
  the wrong start; the missed genes are 2.4%. Wrong starts are split 56%
  "predicted too short" / 44% "too long", mostly within a few codons.
- **The wrong starts are not near-ties.** For *E. coli* and *B. subtilis*, the
  chosen start outscores the annotated start in 90% of the wrong cases, with
  a median total-score gap of 5–7. The RBS component drives most of the gap
  (*B. subtilis*: chosen start's RBS score higher in 62% of cases, median
  +5.0); coding, type and upstream components are roughly neutral. In other
  words, the scoring function is confidently wrong, not undecided.
- **Missed genes are short.** Across the 34 genomes, 9% of the missed
  reference genes are shorter than the 90 bp minimum gene length (impossible
  by design), 26% are 90–149 bp and 44% are 150–299 bp; only 20% are 300 bp
  or longer. Short ORFs carry little coding signal, and lowering the minimum
  would mostly add false positives of the same kind.
- **False positives are mostly short ORFs** (45% under 300 bp); 22% overlap an
  annotated gene on the other strand.

## What was tried

**Re-weighting Prodigal's start scores.** A conditional-logit re-ranker was
trained on the candidate starts of every found gene (features: the coding,
RBS, upstream and type scores and their differences to the best candidate in
the ORF, start-codon type, rank from the 5' end, distance to the most
upstream start, and whether the candidate is Prodigal's choice), with 5-fold
cross-validation by genome. It reproduces Prodigal's accuracy exactly
(90.0% vs 90.1% of found genes) and cannot improve on it: Prodigal's start
score is already an essentially optimal linear combination of its own
features.

**Adding raw sequence context.** The same model with one-hot bases from −30
to +6 around each candidate start (150 extra weights) reaches 90.2%, a gain
of 0.14 points with the same cross-validation: small gains on most genomes,
a 1.2-point loss on *Streptomyces coelicolor*. The signal is real but far
too small to justify a second start model that also breaks Prodigal
equivalence.

## Conclusions and directions

With sequence alone, ab initio start prediction plateaus around 90% on these
references; the remaining disagreements are cases where the sequence signals
(RBS, composition, coding potential) genuinely favour a different start than
the annotation. Closing that gap needs evidence that is not in the genome
being analysed:

- protein homology for candidate starts and short ORFs (how PGAP and
  GeneMarkS-2+ reach their accuracy), or
- experimentally verified start sets (N-terminal proteomics) used to train a
  richer, genome-adaptive start model.

Either is a product-level addition (a protein database or a training corpus
shipped with the tool) rather than a change to the current scoring.
