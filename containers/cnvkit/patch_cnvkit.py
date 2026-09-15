#!/usr/bin/env python3
# Cherry-picks etal/cnvkit PR #1131 (fix for issue #1125: duplicate segment
# boundaries from non-unique pandas index after concat) onto an installed
# cnvkit==0.9.14. Upstream merged 2026-07-22, milestoned for the not-yet-released
# v0.9.15 -- this applies the same two-line fix directly to the installed package.
#
# Run inside the target environment/container:
#   python3 patch_cnvkit.py

import os
import sys

import cnvlib.segmentation as _seg

seg_dir = os.path.dirname(_seg.__file__)
init_path = os.path.join(seg_dir, "__init__.py")
haar_path = os.path.join(seg_dir, "haar.py")


def replace_once(path, old, new, label):
    with open(path) as fh:
        src = fh.read()
    n = src.count(old)
    if n == 0:
        sys.exit(f"FAILED [{label}]: anchor text not found in {path} -- "
                  f"cnvkit source doesn't match what this patch expects.")
    if n > 1:
        sys.exit(f"FAILED [{label}]: anchor text found {n} times in {path}, "
                  f"expected exactly 1 -- refusing to guess which to patch.")
    with open(path, "w") as fh:
        fh.write(src.replace(old, new, 1))
    print(f"OK [{label}]: patched {path}")


# --- __init__.py, hunk 1: CBS/none re-segmentation concat ------------------
replace_once(
    init_path,
    "            segarr = segarr.as_dataframe(pd.concat(newsegs))\n",
    "            # ignore_index: each variants_in_segment result is indexed from 0,\n"
    "            # so concatenating them yields duplicate labels. The baf assignment\n"
    "            # below aligns by label and would broadcast the first segment's BAF\n"
    "            # to every row sharing a label -- a segmentation result must carry a\n"
    "            # unique index (#1125).\n"
    "            segarr = segarr.as_dataframe(pd.concat(newsegs, ignore_index=True))\n",
    "init.py: pd.concat(newsegs, ignore_index=True)",
)

# --- __init__.py, hunk 2: transfer_fields endpoint stretch, loc -> iloc ----
replace_once(
    init_path,
    '    # Avoid chained assignment by directly modifying the underlying DataFrame\n'
    '    segments.data.loc[segments.data.index[0], "start"] = bins_start\n'
    '    segments.data.loc[segments.data.index[-1], "end"] = bins_end\n',
    '    # Stretch the first and last segment endpoints to cover the arm\'s original\n'
    '    # bins. Address by POSITION, not index label: haar concatenates per-arm\n'
    '    # segment tables whose labels can repeat, and `.loc[label]` would broadcast\n'
    '    # the write to every row sharing that label, collapsing distinct segments to\n'
    '    # the full-arm span (#1125).\n'
    '    start_col = segments.data.columns.get_loc("start")\n'
    '    end_col = segments.data.columns.get_loc("end")\n'
    '    segments.data.iloc[0, start_col] = bins_start\n'
    '    segments.data.iloc[-1, end_col] = bins_end\n',
    "init.py: transfer_fields loc -> iloc",
)

# --- haar.py: per-arm segment table concat (not used by this pipeline's -m cbs,
# --- patched anyway for consistency in case haar is ever used) -------------
replace_once(
    haar_path,
    "    segarr = cnarr.as_dataframe(pd.concat(chrom_tables))\n",
    "    # ignore_index: by_arm splits a chromosome into arms, and each arm's segment\n"
    "    # table starts its index at 0, so a plain concat would yield duplicate index\n"
    "    # labels. A segmentation result must carry a unique index for downstream\n"
    "    # label-based operations to behave (#1125).\n"
    "    segarr = cnarr.as_dataframe(pd.concat(chrom_tables, ignore_index=True))\n",
    "haar.py: pd.concat(chrom_tables, ignore_index=True)",
)

print("All three hunks applied successfully.")
