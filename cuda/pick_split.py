#!/usr/bin/env python3
# pick_split.py <adapterremoval.json> <fnr> [budget]
#
# Auto-pick the short/long split length for adna_aligner.sh:
#   reads < L  -> bwa aln / gpualn   (BWA-backtrack; sensitive on short damaged reads)
#   reads >= L -> bwa mem / mem3
#
# Uses AdapterRemoval3's post-processing read-length histogram (output.*.lengths in
# the JSON; free, no FASTQ rescan) plus bwa's own max_diff-vs-length curve to keep the
# `aln` work on the cheap side of the cost cliff. Tree-cost per read ~ C(len,d)*3^d with
# d = bwa_cal_maxdiff(len, 0.02, fnr); cost explodes once d crosses 5.
#
# Rule: largest L such that  (a) all reads < L have max_diff <= 4  (hard cliff guard),
#                            (b) cumulative aln tree-cost(< L) <= budget  (runtime guard).
# Falls back to 50 on any error so the pipeline never breaks.
import json, sys, math
from math import comb

def maxdiff(l, fnr, err=0.02):
    if l < 1: return 0
    el = math.exp(-l*err); s = el; y = 1.0; x = 1
    for k in range(1, 1000):
        y *= l*err; x *= k; s += el*y/x
        if 1.0 - s < fnr: return k
    return 2

def main():
    try:
        js = sys.argv[1]; fnr = float(sys.argv[2])
        budget = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0e14
        d = json.load(open(js))
        # sum every non-discarded output length histogram (covers SE read1, and PE
        # merged/singleton/read1/read2 -- all of which land in the split input)
        out = d.get("output", {})
        hist = {}
        def collect(o):
            if isinstance(o, dict):
                if "lengths" in o and isinstance(o["lengths"], list):
                    for i, c in enumerate(o["lengths"]): hist[i] = hist.get(i, 0) + c
                for k, v in o.items():
                    if k != "lengths": collect(v)
            elif isinstance(o, list):
                for v in o: collect(v)
        for k, v in out.items():
            if k != "discarded": collect(v)
        if not hist: raise ValueError("no length histogram")
        maxlen = max(l for l, c in hist.items() if c)
        cliff = next((l for l in range(1, maxlen+2) if maxdiff(l, fnr) > 4), maxlen+1)  # cliff guard
        def cost(l):
            m = maxdiff(l, fnr); return comb(l, m)*(3**m) if l >= m else 1
        cum = 0; best = 1
        for L in range(1, cliff+1):
            cum += cost(L-1) * hist.get(L-1, 0)   # cost of reads of length L-1 (i.e. < L now)
            if cum <= budget: best = L
            else: break
        print(max(best, 1))
    except Exception as e:
        sys.stderr.write(f"[pick_split] warning: {e}; falling back to 50\n")
        print(50)

if __name__ == "__main__":
    main()
