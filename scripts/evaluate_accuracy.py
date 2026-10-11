#!/usr/bin/env python3
"""Compare Orphos GFF predictions with reference CDS TSVs.

A reference gene counts as found when a prediction shares its 3' end (stop);
it is exact when the 5' end matches too. Predictions without a reference gene
at their 3' end are false positives relative to the annotation.

Usage: evaluate_accuracy.py <ref_dir> <gff_dir> <accession>...
  <ref_dir>/<accession>.cds.tsv and <gff_dir>/<accession>.gff must exist.
"""
import sys, collections


def load_ref(p):
    out = {}
    for line in open(p).read().splitlines()[1:]:
        s, e, st = line.split('\t'); s, e = int(s), int(e)
        key = (e, st) if st == '+' else (s, st)   # 3' end identifies the gene
        out[key] = (s, e, st)
    return out
def load_gff(p):
    out = {}
    for line in open(p):
        if line.startswith('#'): continue
        c = line.rstrip('\n').split('\t')
        if len(c) < 9 or c[2] != 'CDS': continue
        s, e, st = int(c[3]), int(c[4]), c[6]
        key = (e, st) if st == '+' else (s, st)
        attrs = dict(kv.split('=', 1) for kv in c[8].strip(';').split(';') if '=' in kv)
        out[key] = (s, e, st, attrs)
    return out
def evaluate(ref, pred):
    r = collections.Counter(); offsets = collections.Counter(); wrong = []
    for key, (s, e, st) in ref.items():
        if key not in pred:
            r['missed'] += 1; continue
        ps, pe, pst, attrs = pred[key]
        if (ps, pe) == (s, e):
            r['exact'] += 1
        else:
            r['wrong_start'] += 1
            # positive offset = predicted start is downstream (shorter gene)
            off = (ps - s) if st == '+' else (e - pe)
            offsets[off // 3] += 1
            wrong.append((key, s, e, st, ps, pe, attrs))
    r['false_pos'] = sum(1 for k in pred if k not in ref)
    r['ref'] = len(ref); r['pred'] = len(pred)
    return r, offsets, wrong
if __name__ == '__main__':
    tot = collections.Counter()
    print('accession\tref\tpred\tfound(%)\texact(%)\tstart_ok|found(%)\tFP\tmissed\twrong_start')
    for acc in sys.argv[3:]:
        try:
            ref = load_ref(f'{sys.argv[1]}/{acc}.cds.tsv'); pred = load_gff(f'{sys.argv[2]}/{acc}.gff')
        except FileNotFoundError:
            continue
        r, offsets, wrong = evaluate(ref, pred); tot += r
        found = r['exact'] + r['wrong_start']
        print(f"{acc}\t{r['ref']}\t{r['pred']}\t{100*found/r['ref']:.1f}\t{100*r['exact']/r['ref']:.1f}\t{100*r['exact']/max(found,1):.1f}\t{r['false_pos']}\t{r['missed']}\t{r['wrong_start']}")
    found = tot['exact'] + tot['wrong_start']
    print(f"TOTAL\t{tot['ref']}\t{tot['pred']}\t{100*found/tot['ref']:.1f}\t{100*tot['exact']/tot['ref']:.1f}\t{100*tot['exact']/max(found,1):.1f}\t{tot['false_pos']}\t{tot['missed']}\t{tot['wrong_start']}")
