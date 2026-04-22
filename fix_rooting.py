#!/usr/bin/env python3
"""
Re-roots tree.nwk files that Biopython couldn't parse (apostrophes in names).
Uses dendropy which handles FastTree output robustly.
Only processes trees where the reference is not yet at the root.
"""

import os
import re
import dendropy
from Bio import SeqIO

BASE_DIR   = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "split_output")


def sanitize_name(name: str) -> str:
    for ch in "(),;:[]'\"`\t ":
        name = name.replace(ch, "_")
    return name


def get_ref_id(ref_path: str) -> str:
    return sanitize_name(next(SeqIO.parse(ref_path, "fasta")).id)


fixed = skipped = failed = 0

for subtype in sorted(os.listdir(OUTPUT_DIR)):
    sdir = os.path.join(OUTPUT_DIR, subtype)
    if not os.path.isdir(sdir): continue
    for segment in sorted(os.listdir(sdir)):
        segdir = os.path.join(sdir, segment)
        if not os.path.isdir(segdir): continue
        for epoch in sorted(os.listdir(segdir)):
            epoch_dir = os.path.join(segdir, epoch)
            tree_path = os.path.join(epoch_dir, "tree.nwk")
            ref_path  = os.path.join(epoch_dir, "reference.fasta")

            if not os.path.exists(tree_path) or not os.path.exists(ref_path):
                continue

            ref_id = get_ref_id(ref_path)

            # Check if already rooted correctly using simple string check
            with open(tree_path) as f:
                nwk = f.read()

            # Sanitize the raw NWK string (replace apostrophes etc.)
            sanitized_nwk = nwk
            for ch in "'`":
                sanitized_nwk = sanitized_nwk.replace(ch, "_")

            # Check if ref is already a direct child of root
            # (appears right after opening paren or comma at depth 0)
            ref_at_root = bool(re.match(
                r'\(' + re.escape(ref_id.split("|")[0].replace("'", "_")),
                sanitized_nwk
            ))

            if ref_at_root:
                skipped += 1
                continue

            try:
                # Write sanitized NWK to temp file for dendropy
                sanitized_path = tree_path + ".tmp"
                with open(sanitized_path, "w") as f:
                    f.write(sanitized_nwk)

                tree = dendropy.Tree.get(path=sanitized_path, schema="newick",
                                         preserve_underscores=True)
                os.remove(sanitized_path)

                # Find reference leaf
                ref_node = None
                for leaf in tree.leaf_node_iter():
                    if leaf.taxon and leaf.taxon.label and \
                       (leaf.taxon.label == ref_id or
                        leaf.taxon.label.startswith(ref_id.split("|")[0])):
                        ref_node = leaf
                        break

                if ref_node is None:
                    failed += 1
                    print(f"  [REF-NOT-FOUND] {subtype}/{segment}/{epoch}")
                    continue

                tree.reroot_at_edge(ref_node.edge, update_bipartitions=False)
                tree.write(path=tree_path, schema="newick")
                fixed += 1
                print(f"  [FIXED] {subtype}/{segment}/{epoch}")

            except Exception as e:
                failed += 1
                if os.path.exists(sanitized_path):
                    os.remove(sanitized_path)
                print(f"  [FAIL]  {subtype}/{segment}/{epoch}  — {e}")

print(f"\nFixed  : {fixed}")
print(f"Already rooted : {skipped}")
print(f"Failed : {failed}")
