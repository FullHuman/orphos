# Prodigal reference outputs

Inputs and the GFF output of Prodigal for them, compared with the output of
Orphos line by line in `tests/prodigal_reference_outputs.rs`.

The expected files were written by a build of Prodigal's master branch
(commit c1e2d36), not by the 2.6.3 release: the release reads uninitialised
memory in its Shine-Dalgarno searches and applies its short-gene penalty to
every metagenomic gene, both fixed upstream since.

| Input | What it is |
|---|---|
| `contigs.fasta` | 90 pieces of 60 to 9000 bases cut from seven genomes, half of them reverse complemented, as in a metagenome. Some are too short to hold a gene. |
| `gaps.fasta` | Scaffolds with runs of `N` of several lengths, inside and at both ends; lower case; ambiguity codes; a sequence of `N` only; a 75 base sequence; a contig whose last gene starts 14 bases from its end; and sequences of two bases and of one. |
| `draft.fasta` | 220 kb of *H. pylori* 26695 with two gaps, long enough to train on. |

Each expected file is `prodigal -q -i <input> -f gff -o <expected>` with the
options the test lists for it. To check the expected files against a
Prodigal binary, run the tests with `ORPHOS_PRODIGAL_BIN` set to it.
