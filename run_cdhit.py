#!/usr/bin/env python3
"""
Runs cd-hit-est on sequences.fasta in every epoch folder under split_output/.
Output per epoch: clustered.fasta + clustered.fasta.clstr
Identity threshold: 100% (removes exact duplicates only)
"""

import os
import subprocess

BASE_DIR   = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "split_output")

total = skipped = succeeded = failed = 0

for subtype in sorted(os.listdir(OUTPUT_DIR)):
    sdir = os.path.join(OUTPUT_DIR, subtype)
    if not os.path.isdir(sdir): continue

    for segment in sorted(os.listdir(sdir)):
        segdir = os.path.join(sdir, segment)
        if not os.path.isdir(segdir): continue

        print(f"\n[{subtype}] {segment}")

        for epoch in sorted(os.listdir(segdir)):
            epoch_dir = os.path.join(segdir, epoch)
            if not os.path.isdir(epoch_dir): continue

            seq_in = os.path.join(epoch_dir, "sequences.fasta")
            clust  = os.path.join(epoch_dir, "clustered.fasta")

            if not os.path.exists(seq_in):
                continue

            total += 1

            if os.path.exists(clust):
                skipped += 1
                print(f"  [--] {epoch}  (already done)")
                continue

            n_in = sum(1 for l in open(seq_in) if l.startswith(">"))

            try:
                subprocess.run(
                    ["cd-hit-est",
                     "-i", seq_in,
                     "-o", clust,
                     "-c", "1.0",
                     "-n", "8",
                     "-M", "0",
                     "-T", "0",
                     "-d", "0"],
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )
                n_out = sum(1 for l in open(clust) if l.startswith(">"))
                removed = n_in - n_out
                succeeded += 1
                print(f"  [OK] {epoch}  {n_in:,} → {n_out:,} seqs  ({removed:,} duplicates removed)")

            except subprocess.CalledProcessError as e:
                failed += 1
                print(f"  [FAIL] {epoch}  — {e}")

print(f"\n{'='*60}")
print(f" Done.  {succeeded} clustered | {skipped} skipped | {failed} failed | {total} total")
print(f"{'='*60}")
