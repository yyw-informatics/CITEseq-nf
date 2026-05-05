#!/usr/bin/env python
"""Select the winning decontamination branch and promote it to canonical paths.

Copies the branch-suffixed h5ad files (e.g. adata_qc_decontx.h5ad) to the
unsuffixed canonical paths (adata_qc.h5ad) that all downstream steps consume.
The non-selected branch files are left in place for reference.

Usage:
    python scripts/admin/select_branch.py decontx
    python scripts/admin/select_branch.py soupx
    python scripts/admin/select_branch.py decontx --dry-run
"""
import sys
import os
import shutil
import argparse
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from config import PROCESSED_DIR, RESULTS_PHASE1_QC, PROJECT_ROOT

# Files to promote: (suffix_template, canonical_path)
# suffix_template uses {branch} placeholder
BRANCH_FILES = [
    (os.path.join(PROCESSED_DIR, "adata_qc_{branch}.h5ad"),
     os.path.join(PROCESSED_DIR, "adata_qc.h5ad")),
    (os.path.join(PROCESSED_DIR, "adata_adt_raw_{branch}.h5ad"),
     os.path.join(PROCESSED_DIR, "adata_adt_raw.h5ad")),
]

BRANCH_RESULTS = [
    (os.path.join(RESULTS_PHASE1_QC, "doublet_comparison_{branch}.csv"),
     os.path.join(RESULTS_PHASE1_QC, "doublet_comparison.csv")),
]


def select_branch(branch, dry_run=False):
    print(f"=== Selecting branch: {branch} ===")

    all_files = BRANCH_FILES + BRANCH_RESULTS
    # Verify all source files exist before doing anything
    for src_template, _ in all_files:
        src = src_template.format(branch=branch)
        if not os.path.exists(src):
            print(f"  ERROR: source not found: {src}")
            sys.exit(1)

    for src_template, dst in all_files:
        src = src_template.format(branch=branch)
        src_size = os.path.getsize(src) / (1024 * 1024)

        if os.path.exists(dst):
            dst_size = os.path.getsize(dst) / (1024 * 1024)
            action = f"overwrite ({dst_size:.0f} MB -> {src_size:.0f} MB)"
        else:
            action = f"create ({src_size:.0f} MB)"

        rel_src = os.path.relpath(src, PROJECT_ROOT)
        rel_dst = os.path.relpath(dst, PROJECT_ROOT)
        print(f"  {rel_src} -> {rel_dst}  [{action}]")

        if not dry_run:
            shutil.copy2(src, dst)

    # Write a breadcrumb so we know which branch is active
    breadcrumb = os.path.join(PROCESSED_DIR, "selected_branch.txt")
    if not dry_run:
        with open(breadcrumb, "w") as f:
            f.write(f"{branch}\n")
            f.write(f"Selected: {datetime.now().isoformat()}\n")

    if dry_run:
        print("\n  (dry run — no files copied)")
    else:
        print(f"\n  Done. Canonical files now point to '{branch}' branch.")
        print(f"  Breadcrumb: {os.path.relpath(breadcrumb, PROJECT_ROOT)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Promote a decontamination branch to canonical paths")
    parser.add_argument("branch", choices=["decontx", "soupx"],
                        help="Branch to select as canonical")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be copied without doing it")
    args = parser.parse_args()
    select_branch(args.branch, args.dry_run)
