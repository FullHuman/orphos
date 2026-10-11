#!/usr/bin/env python3
"""Convert an NCBI feature table into the start/end/strand TSV used by the
reference-annotation tests.

Fetch a feature table with:
  curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=ft&retmode=text"

Only complete, single-interval, non-pseudo CDS features are kept, so joined
(frameshifted) and pseudogene records do not count as genes to predict.

Usage: refseq_cds_tsv.py <accession.ft> <accession.cds.tsv>
"""
import sys
def parse(path):
    feats = []  # (intervals, qualifiers)
    cur = None
    for line in open(path):
        line = line.rstrip('\n')
        if not line or line.startswith('>'):
            continue
        parts = line.split('\t')
        if len(parts) >= 3 and parts[2]:
            cur = ([(parts[0], parts[1])], parts[2], set())
            feats.append(cur)
        elif len(parts) >= 2 and parts[0] and cur is not None and not parts[0].startswith('\t'):
            cur[0].append((parts[0], parts[1]))
        elif len(parts) >= 4 and cur is not None:
            cur[2].add(parts[3])
    rows = []
    for intervals, kind, quals in feats:
        if kind != 'CDS' or 'pseudo' in quals or len(intervals) != 1:
            continue
        a, b = intervals[0]
        if '<' in a or '>' in a or '<' in b or '>' in b:
            continue  # partial
        a, b = int(a), int(b)
        if a <= b:
            rows.append((a, b, '+'))
        else:
            rows.append((b, a, '-'))
    return rows
rows = parse(sys.argv[1])
with open(sys.argv[2], 'w') as f:
    f.write('start\tend\tstrand\n')
    for r in rows:
        f.write(f'{r[0]}\t{r[1]}\t{r[2]}\n')
print(len(rows))
