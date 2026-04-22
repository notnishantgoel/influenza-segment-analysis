#!/usr/bin/env python3
"""
Builds a phylogenetic tree for every epoch folder that has a
nextclade.aligned.fasta file.

Process per epoch:
  1. Run FastTree (GTR + CAT model) on nextclade.aligned.fasta
  2. Re-root the unrooted tree at the epoch reference sequence
     so all branches show divergence from that epoch's vaccine strain
  3. Write the rooted tree to tree.nwk

FastTree flags used:
  -nt          nucleotide sequences
  -gtr         GTR substitution model (standard for flu evolution)
  -fastest     fastest heuristic (good for large epochs)
  -nosupport   skip bootstrap (saves time; add later if needed)
  -quiet       suppress progress output

Idempotent: skips epochs where tree.nwk already exists.
Requires: FastTree (brew install brewsci/bio/fasttree)
          biopython (pip install biopython)
"""

import os
import subprocess
import time
from Bio import Phylo, SeqIO

BASE_DIR   = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "split_output")
FASTTREE   = "FastTree"

total = skipped = succeeded = failed = no_aln = 0


def get_ref_name(reference_fasta: str) -> str:
    """Return the sanitized sequence ID of the epoch reference."""
    rec = next(SeqIO.parse(reference_fasta, "fasta"))
    return sanitize_name(rec.id)


def sanitize_name(name: str) -> str:
    """Replace Newick special characters with underscores."""
    for ch in "(),;:[]'\"`\t ":
        name = name.replace(ch, "_")
    return name


def dedup_fasta(fasta_path: str, out_path: str) -> int:
    """
    Write fasta_path to out_path with:
      - Newick special characters removed from sequence names
      - Duplicate names removed (first occurrence kept)
    Returns number of duplicates removed.
    """
    seen = set()
    removed = 0
    with open(out_path, "w") as out:
        for rec in SeqIO.parse(fasta_path, "fasta"):
            rec.id = sanitize_name(rec.id)
            rec.description = rec.id
            if rec.id in seen:
                removed += 1
                continue
            seen.add(rec.id)
            SeqIO.write([rec], out, "fasta")
    return removed


def root_at(tree_path: str, ref_id: str) -> str:
    """
    Re-root the Newick file at the leaf matching ref_id, write back.
    Returns 'rooted', 'ref-not-found', or 'parse-error' (file kept as-is).
    """
    try:
        tree = Phylo.read(tree_path, "newick")
    except Exception:
        return "parse-error"
    targets = [c for c in tree.find_clades(terminal=True)
               if c.name and (c.name == ref_id or c.name.startswith(ref_id))]
    if not targets:
        return "ref-not-found"
    tree.root_with_outgroup(targets[0])
    Phylo.write(tree, tree_path, "newick")
    return "rooted"


for subtype in sorted(os.listdir(OUTPUT_DIR)):
    sdir = os.path.join(OUTPUT_DIR, subtype)
    if not os.path.isdir(sdir):
        continue

    for segment in sorted(os.listdir(sdir)):
        segdir = os.path.join(sdir, segment)
        if not os.path.isdir(segdir):
            continue

        print(f"\n[{subtype}] {segment}")

        for epoch in sorted(os.listdir(segdir)):
            epoch_dir = os.path.join(segdir, epoch)
            if not os.path.isdir(epoch_dir):
                continue

            aligned   = os.path.join(epoch_dir, "nextclade.aligned.fasta")
            reference = os.path.join(epoch_dir, "reference.fasta")
            tree_out  = os.path.join(epoch_dir, "tree.nwk")

            if not os.path.exists(aligned):
                no_aln += 1
                continue

            total += 1

            if os.path.exists(tree_out):
                skipped += 1
                n = sum(1 for l in open(aligned) if l.startswith(">"))
                print(f"  [--]   {epoch}  (already done)")
                continue

            n_seqs = sum(1 for l in open(aligned) if l.startswith(">"))

            # Deduplicate alignment (reference may appear in both ref + clustered)
            dedup_path = aligned + ".dedup.fasta"
            removed = dedup_fasta(aligned, dedup_path)
            input_fasta = dedup_path

            t0 = time.time()
            try:
                result = subprocess.run(
                    [FASTTREE, "-nt", "-gtr", "-fastest", "-nosupport", "-quiet",
                     input_fasta],
                    check=True,
                    capture_output=True,
                    text=True,
                )
                nwk_string = result.stdout
                elapsed = time.time() - t0

                # Write raw NWK first, then re-root at epoch reference
                with open(tree_out, "w") as f:
                    f.write(nwk_string)

                if os.path.exists(reference):
                    ref_id = get_ref_name(reference)
                    rooted = root_at(tree_out, ref_id)
                else:
                    rooted = "unrooted"

                succeeded += 1
                dedup_note = f"  (-{removed} dupes)" if removed else ""
                print(f"  [OK]   {epoch}  {n_seqs:,} seqs  [{rooted}]{dedup_note}  ({elapsed:.0f}s)")

            except subprocess.CalledProcessError as e:
                elapsed = time.time() - t0
                failed += 1
                if os.path.exists(tree_out):
                    os.remove(tree_out)
                print(f"  [FAIL] {epoch}  — {e.stderr[:80] if e.stderr else e}  ({elapsed:.0f}s)")
            finally:
                if os.path.exists(dedup_path):
                    os.remove(dedup_path)

print(f"\n{'='*65}")
print(f" Done.")
print(f"   Trees built : {succeeded}")
print(f"   Skipped     : {skipped}  (already done)")
print(f"   No alignment: {no_aln}  (New Caledonia — excluded)")
print(f"   Failed      : {failed}")
print(f"   Total       : {total}")
print(f"{'='*65}")
