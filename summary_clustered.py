#!/usr/bin/env python3
"""
Summary of clustered sequences per epoch for H1N1 and H3N2.
Reports sequence counts before (sequences.fasta) and after (clustered.fasta)
CD-HIT deduplication, and the reduction percentage.
"""

import os
from collections import defaultdict

BASE_DIR   = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "split_output")


def count_seqs(path):
    if not os.path.exists(path):
        return 0
    return sum(1 for l in open(path) if l.startswith(">"))


for subtype in sorted(os.listdir(OUTPUT_DIR)):
    sdir = os.path.join(OUTPUT_DIR, subtype)
    if not os.path.isdir(sdir):
        continue

    print(f"\n{'='*75}")
    print(f"  {subtype}")
    print(f"{'='*75}")

    subtype_orig = subtype_clust = 0

    for segment in sorted(os.listdir(sdir)):
        segdir = os.path.join(sdir, segment)
        if not os.path.isdir(segdir):
            continue

        print(f"\n  [{segment}]")
        print(f"  {'Epoch':<45} {'Original':>9} {'Clustered':>10} {'Removed':>8} {'Reduction':>10}")
        print(f"  {'─'*45} {'─'*9} {'─'*10} {'─'*8} {'─'*10}")

        seg_orig = seg_clust = 0

        for epoch in sorted(os.listdir(segdir)):
            epoch_dir = os.path.join(segdir, epoch)
            if not os.path.isdir(epoch_dir):
                continue

            orig  = count_seqs(os.path.join(epoch_dir, "sequences.fasta"))
            clust = count_seqs(os.path.join(epoch_dir, "clustered.fasta"))
            removed = orig - clust
            pct = (removed / orig * 100) if orig > 0 else 0.0

            print(f"  {epoch:<45} {orig:>9,} {clust:>10,} {removed:>8,} {pct:>9.1f}%")

            seg_orig  += orig
            seg_clust += clust

        seg_removed = seg_orig - seg_clust
        seg_pct = (seg_removed / seg_orig * 100) if seg_orig > 0 else 0.0
        print(f"  {'─'*45} {'─'*9} {'─'*10} {'─'*8} {'─'*10}")
        print(f"  {'TOTAL ' + segment:<45} {seg_orig:>9,} {seg_clust:>10,} {seg_removed:>8,} {seg_pct:>9.1f}%")

        subtype_orig  += seg_orig
        subtype_clust += seg_clust

    sub_removed = subtype_orig - subtype_clust
    sub_pct = (sub_removed / subtype_orig * 100) if subtype_orig > 0 else 0.0
    print(f"\n  {'═'*75}")
    print(f"  {'TOTAL ' + subtype:<45} {subtype_orig:>9,} {subtype_clust:>10,} {sub_removed:>8,} {sub_pct:>9.1f}%")

print(f"\n{'='*75}\n Done.\n{'='*75}")
