#!/usr/bin/env python3
"""
ICassigner.py v0.1.0

Conservative phylogeny-guided assignment of Acinetobacter baumannii
International Clones (ICs) using a core genome tree and a metadata CSV.

Features
- Conservative tree propagation of IC labels (default: UA -> IC only if support>=10 and prop>=0.90)
- Optional plots:
    1) IC counts before vs after
    2) support_n histogram (inferred only)
    3) support_prop histogram (inferred only)
    4) confusion matrices as heatmaps (ST-expected IC vs assigned IC; group columns vs IC)
- Cross-checks:
    * Pasteur ST (ST1, ST2...) -> expected IC
    * Oxford ST (common anchors) -> expected IC
    * Any user-supplied grouping column (cgMLST group, BAPS, etc.) -> confusion vs IC

Dependencies
- ete3
- matplotlib (optional; only needed if --plots)
- numpy (optional; only needed if --plots)

Example
python ICassigner.py --tree Acinetobacter-coreML.nwk \
  --metadata metadata.csv \
  --tip_col sample_id --ic_col IC \
  --plots --outdir outputs_icassigner \
  --pasteur_st_col ST_Pasteur --oxford_st_col ST_Oxford \
  --group_cols cgMLST_group,BAPS
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from ete3 import Tree


# --- Canonical ST → IC anchors (sanity check, not primary lineage definition) ---
PASTEUR_ST_TO_IC = {
    "ST1": "IC1",
    "ST2": "IC2",
    "ST3": "IC3",
    "ST15": "IC4",
    "ST79": "IC5",
    "ST78": "IC6",
    "ST25": "IC7",
    "ST10": "IC8",
    "ST85": "IC9",
}

# Oxford ST mapping is less one-to-one; these are common anchors used as checks
OXFORD_ST_TO_IC = {
    "ST231": "IC1",
    "ST208": "IC2",
    "ST451": "IC2",
    "ST348": "IC3",
    "ST15": "IC4",
    "ST79": "IC5",
    "ST78": "IC6",
    "ST25": "IC7",
    "ST85": "IC9",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Conservative IC assignment from a phylogenetic tree (pandas-free).")
    p.add_argument("--tree", required=True, help="Newick tree file (e.g. RAxML core genome tree)")
    p.add_argument("--metadata", required=True, help="Metadata CSV containing sample IDs and IC labels")
    p.add_argument("--tip_col", default="sample_id", help="Column for tree tip/sample IDs")
    p.add_argument("--ic_col", default="IC", help="Column for International Clone labels")
    p.add_argument("--ua_label", default="UA", help='Label used for "unassigned" (default: UA)')

    p.add_argument("--min_support", type=int, default=10, help="Minimum labelled neighbours required (default: 10)")
    p.add_argument("--min_prop", type=float, default=0.90, help="Minimum majority proportion required (default: 0.90)")

    p.add_argument("--output", default="metadata_with_conservative_IC.csv", help="Output CSV filename")
    p.add_argument("--outdir", default="outputs_icassigner", help="Output directory (default: outputs_icassigner)")
    p.add_argument("--plots", action="store_true", help="Generate plots (requires matplotlib + numpy)")

    p.add_argument("--pasteur_st_col", default=None, help="Column containing Pasteur MLST ST (e.g., ST_Pasteur)")
    p.add_argument("--oxford_st_col", default=None, help="Column containing Oxford MLST ST (e.g., ST_Oxford)")
    p.add_argument("--group_cols", default=None,
                   help="Comma-separated list of additional group columns (e.g., cgMLST_group,BAPS)")
    p.add_argument("--max_groups_plot", type=int, default=30,
                   help="Max categories to show on confusion plots for group_cols (default: 30). Full CSV is always written.")

    p.add_argument("--force_overwrite", action="store_true",
                   help="Overwrite existing non-UA IC labels with tree-inferred labels (NOT recommended).")

    p.add_argument("--loglevel", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                   help="Logging level (default: INFO)")
    return p.parse_args()


def setup_logging(level: str) -> None:
    logging.basicConfig(level=getattr(logging, level), format="%(asctime)s | %(levelname)s | %(message)s")


def safe_makedirs(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def read_csv_rows(path: str) -> Tuple[List[Dict[str, str]], List[str]]:
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        rows = [dict(r) for r in reader]
        return rows, reader.fieldnames or []


def write_csv_rows(path: str, rows: List[Dict[str, str]], fieldnames: List[str]) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def normalize_st(x: Optional[str]) -> Optional[str]:
    if x is None:
        return None
    s = str(x).strip()
    if not s or s.lower() == "na":
        return None
    if s.upper().startswith("ST"):
        return "ST" + s[2:].strip()
    if s.isdigit():
        return "ST" + s
    return s


@dataclass
class InferenceResult:
    ic: str
    support_n: int
    majority_n: int
    majority_prop: float
    node_size: int


def infer_ic_from_tree(
    tree: Tree,
    tip: str,
    ic_map: Dict[str, Optional[str]],
    ua_label: str,
    min_support: int,
    min_prop: float,
) -> InferenceResult:
    node = tree & tip
    while node is not None:
        tips = node.get_leaf_names()
        labels = [ic_map.get(t) for t in tips if ic_map.get(t) not in [None, ua_label, ""]]
        if len(labels) >= min_support:
            counts = Counter(labels)
            ic, majority_n = counts.most_common(1)[0]
            total = sum(counts.values())
            prop = majority_n / total if total else 0.0
            if prop >= min_prop:
                return InferenceResult(ic=ic, support_n=total, majority_n=majority_n, majority_prop=prop, node_size=len(tips))
        node = node.up
    return InferenceResult(ic=ua_label, support_n=0, majority_n=0, majority_prop=0.0, node_size=0)


def value_counts(values: List[str]) -> Counter:
    c = Counter()
    for v in values:
        c[v] += 1
    return c


def write_summary(rows: List[Dict[str, str]], ic_col: str, out_ic_col: str, ua_label: str, outdir: str) -> None:
    safe_makedirs(outdir)
    before = [r.get(ic_col, "") or "" for r in rows]
    after = [r.get(out_ic_col, "") or "" for r in rows]

    total = len(rows)
    ua_before = sum(1 for v in before if v == ua_label)
    ua_after = sum(1 for v in after if v == ua_label)
    inferred = sum(1 for b, a in zip(before, after) if b == ua_label and a != ua_label)

    before_counts = value_counts(before)
    after_counts = value_counts(after)

    lines = [
        f"Total isolates: {total}",
        f"UA before: {ua_before}",
        f"UA after: {ua_after}",
        f"Inferred (UA -> IC): {inferred}",
        "",
        "Counts before (top 20):",
        *[f"  {k}: {v}" for k, v in before_counts.most_common(20)],
        "",
        "Counts after (top 20):",
        *[f"  {k}: {v}" for k, v in after_counts.most_common(20)],
        "",
    ]
    with open(os.path.join(outdir, "ICassigner_summary.txt"), "w") as f:
        f.write("\n".join(lines))


def crosstab(rows: List[Dict[str, str]], row_key: str, col_key: str, mask_fn=None) -> Tuple[List[str], List[str], Dict[Tuple[str, str], int]]:
    """Return (row_levels, col_levels, counts[(r,c)])."""
    counts = defaultdict(int)
    row_levels_set, col_levels_set = set(), set()
    for r in rows:
        if mask_fn and not mask_fn(r):
            continue
        rv = (r.get(row_key, "") or "").strip()
        cv = (r.get(col_key, "") or "").strip()
        if rv == "" or cv == "":
            continue
        counts[(rv, cv)] += 1
        row_levels_set.add(rv)
        col_levels_set.add(cv)
    row_levels = sorted(row_levels_set)
    col_levels = sorted(col_levels_set)
    return row_levels, col_levels, counts


def write_crosstab_csv(out_path: str, row_levels: List[str], col_levels: List[str], counts: Dict[Tuple[str, str], int],
                       row_label: str, col_label: str) -> None:
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([row_label] + col_levels)
        for rv in row_levels:
            w.writerow([rv] + [counts.get((rv, cv), 0) for cv in col_levels])


def plot_confusion_heatmap(out_png: str, out_pdf: str, title: str,
                           row_levels: List[str], col_levels: List[str], counts: Dict[Tuple[str, str], int],
                           max_rows: int, max_cols: int) -> None:
    import numpy as np
    import matplotlib.pyplot as plt

    # compute marginals for truncation
    row_totals = {rv: sum(counts.get((rv, cv), 0) for cv in col_levels) for rv in row_levels}
    col_totals = {cv: sum(counts.get((rv, cv), 0) for rv in row_levels) for cv in col_levels}

    if len(row_levels) > max_rows:
        row_levels = [rv for rv, _ in sorted(row_totals.items(), key=lambda x: x[1], reverse=True)[:max_rows]]
    if len(col_levels) > max_cols:
        col_levels = [cv for cv, _ in sorted(col_totals.items(), key=lambda x: x[1], reverse=True)[:max_cols]]

    data = np.array([[counts.get((rv, cv), 0) for cv in col_levels] for rv in row_levels], dtype=float)

    plt.figure(figsize=(max(6, 0.3 * len(col_levels)), max(4, 0.3 * len(row_levels))))
    plt.imshow(data, aspect="auto")
    plt.title(title)
    plt.xticks(range(len(col_levels)), col_levels, rotation=90)
    plt.yticks(range(len(row_levels)), row_levels)
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()


def plot_ic_before_after(rows: List[Dict[str, str]], ic_col: str, out_ic_col: str, outdir: str) -> None:
    import matplotlib.pyplot as plt

    before = [(r.get(ic_col, "") or "NA").strip() for r in rows]
    after = [(r.get(out_ic_col, "") or "NA").strip() for r in rows]
    all_ics = sorted(set(before) | set(after))

    before_counts = Counter(before)
    after_counts = Counter(after)

    x = list(range(len(all_ics)))
    bvals = [before_counts.get(ic, 0) for ic in all_ics]
    avals = [after_counts.get(ic, 0) for ic in all_ics]

    plt.figure()
    plt.bar([i - 0.2 for i in x], bvals, width=0.4, label="Before")
    plt.bar([i + 0.2 for i in x], avals, width=0.4, label="After")
    plt.xticks(x, all_ics, rotation=90)
    plt.ylabel("Isolates")
    plt.title("IC counts before vs after conservative tree inference")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "IC_counts_before_after.png"), dpi=300)
    plt.savefig(os.path.join(outdir, "IC_counts_before_after.pdf"))
    plt.close()


def plot_histograms(rows: List[Dict[str, str]], ic_col: str, out_ic_col: str, ua_label: str, outdir: str) -> None:
    import matplotlib.pyplot as plt

    inferred_support_n = []
    inferred_support_prop = []

    for r in rows:
        before = (r.get(ic_col, "") or "").strip()
        after = (r.get(out_ic_col, "") or "").strip()
        if before == ua_label and after != ua_label:
            try:
                inferred_support_n.append(float(r.get("IC_tree_support_n", "0") or 0))
            except Exception:
                pass
            try:
                inferred_support_prop.append(float(r.get("IC_tree_support_prop", "0") or 0))
            except Exception:
                pass

    if inferred_support_n:
        plt.figure()
        plt.hist(inferred_support_n, bins=30)
        plt.xlabel("Labelled neighbours used (support_n)")
        plt.ylabel("Count")
        plt.title("Support used for inferred IC calls")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "IC_inference_support_hist.png"), dpi=300)
        plt.savefig(os.path.join(outdir, "IC_inference_support_hist.pdf"))
        plt.close()

    if inferred_support_prop:
        plt.figure()
        plt.hist(inferred_support_prop, bins=30)
        plt.xlabel("Majority support proportion (support_prop)")
        plt.ylabel("Count")
        plt.title("Majority support for inferred IC calls")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "IC_inference_prop_hist.png"), dpi=300)
        plt.savefig(os.path.join(outdir, "IC_inference_prop_hist.pdf"))
        plt.close()


def main() -> None:
    args = parse_args()
    setup_logging(args.loglevel)
    safe_makedirs(args.outdir)

    logging.info("Reading tree: %s", args.tree)
    tree = Tree(args.tree, format=1)
    tree_tips = set(tree.get_leaf_names())

    logging.info("Reading metadata: %s", args.metadata)
    rows, fieldnames = read_csv_rows(args.metadata)

    if args.tip_col not in fieldnames:
        raise ValueError(f"--tip_col '{args.tip_col}' not found in metadata columns.")
    if args.ic_col not in fieldnames:
        raise ValueError(f"--ic_col '{args.ic_col}' not found in metadata columns.")

    # IC map keyed by tip/sample id
    ic_map: Dict[str, Optional[str]] = {}
    for r in rows:
        tip = (r.get(args.tip_col, "") or "").strip()
        ic = (r.get(args.ic_col, "") or "").strip()
        ic_map[tip] = ic if ic else None

    # Identify metadata tips missing from tree
    missing_in_tree = sorted([tip for tip in ic_map.keys() if tip and tip not in tree_tips])
    if missing_in_tree:
        logging.warning("Found %d metadata tips not present in tree. These cannot be inferred.", len(missing_in_tree))
        with open(os.path.join(args.outdir, "tips_missing_in_tree.txt"), "w") as f:
            f.write("\n".join(missing_in_tree) + "\n")

    out_ic_col = "IC_tree_conservative"
    # add columns (ensure they are in output fieldnames)
    extra_cols = [
        out_ic_col,
        "IC_tree_support_n",
        "IC_tree_majority_n",
        "IC_tree_support_prop",
        "IC_tree_node_size",
    ]

    for r in rows:
        tip = (r.get(args.tip_col, "") or "").strip()
        original_ic = (r.get(args.ic_col, "") or "").strip()

        if (not args.force_overwrite) and original_ic not in ["", args.ua_label]:
            # Keep known label, leave inference stats blank
            r[out_ic_col] = original_ic
            r["IC_tree_support_n"] = ""
            r["IC_tree_majority_n"] = ""
            r["IC_tree_support_prop"] = ""
            r["IC_tree_node_size"] = ""
            continue

        if tip not in tree_tips:
            # Can't infer
            r[out_ic_col] = original_ic if original_ic else args.ua_label
            r["IC_tree_support_n"] = "0"
            r["IC_tree_majority_n"] = "0"
            r["IC_tree_support_prop"] = "0.0"
            r["IC_tree_node_size"] = "0"
            continue

        res = infer_ic_from_tree(
            tree=tree,
            tip=tip,
            ic_map=ic_map,
            ua_label=args.ua_label,
            min_support=args.min_support,
            min_prop=args.min_prop,
        )
        r[out_ic_col] = res.ic
        r["IC_tree_support_n"] = str(res.support_n)
        r["IC_tree_majority_n"] = str(res.majority_n)
        r["IC_tree_support_prop"] = f"{res.majority_prop:.6f}"
        r["IC_tree_node_size"] = str(res.node_size)

    # Cross-checks: expected IC from Pasteur/Oxford ST
    if args.pasteur_st_col and args.pasteur_st_col in fieldnames:
        for r in rows:
            st = normalize_st(r.get(args.pasteur_st_col))
            exp = PASTEUR_ST_TO_IC.get(st) if st else None
            r["IC_expected_from_PasteurST"] = exp or ""
            assigned = (r.get(out_ic_col, "") or "").strip()
            r["IC_PasteurST_conflict"] = "1" if (exp and assigned and assigned != args.ua_label and exp != assigned) else "0"
        extra_cols += ["IC_expected_from_PasteurST", "IC_PasteurST_conflict"]

        mask_fn = lambda rr: (rr.get("IC_expected_from_PasteurST", "") or "").strip() != "" and (rr.get(out_ic_col, "") or "").strip() != args.ua_label
        rl, cl, ct = crosstab(rows, "IC_expected_from_PasteurST", out_ic_col, mask_fn=mask_fn)
        out_csv = os.path.join(args.outdir, "confusion_expectedIC_PasteurST_vs_assignedIC.csv")
        write_crosstab_csv(out_csv, rl, cl, ct, "Expected_IC (Pasteur ST)", "Assigned_IC (tree)")

        if args.plots:
            plot_confusion_heatmap(
                out_png=os.path.join(args.outdir, "confusion_expectedIC_PasteurST_vs_assignedIC.png"),
                out_pdf=os.path.join(args.outdir, "confusion_expectedIC_PasteurST_vs_assignedIC.pdf"),
                title="Pasteur ST-expected IC vs Tree-assigned IC",
                row_levels=rl, col_levels=cl, counts=ct,
                max_rows=50, max_cols=50,
            )

    if args.oxford_st_col and args.oxford_st_col in fieldnames:
        for r in rows:
            st = normalize_st(r.get(args.oxford_st_col))
            exp = OXFORD_ST_TO_IC.get(st) if st else None
            r["IC_expected_from_OxfordST"] = exp or ""
            assigned = (r.get(out_ic_col, "") or "").strip()
            r["IC_OxfordST_conflict"] = "1" if (exp and assigned and assigned != args.ua_label and exp != assigned) else "0"
        extra_cols += ["IC_expected_from_OxfordST", "IC_OxfordST_conflict"]

        mask_fn = lambda rr: (rr.get("IC_expected_from_OxfordST", "") or "").strip() != "" and (rr.get(out_ic_col, "") or "").strip() != args.ua_label
        rl, cl, ct = crosstab(rows, "IC_expected_from_OxfordST", out_ic_col, mask_fn=mask_fn)
        out_csv = os.path.join(args.outdir, "confusion_expectedIC_OxfordST_vs_assignedIC.csv")
        write_crosstab_csv(out_csv, rl, cl, ct, "Expected_IC (Oxford ST)", "Assigned_IC (tree)")

        if args.plots:
            plot_confusion_heatmap(
                out_png=os.path.join(args.outdir, "confusion_expectedIC_OxfordST_vs_assignedIC.png"),
                out_pdf=os.path.join(args.outdir, "confusion_expectedIC_OxfordST_vs_assignedIC.pdf"),
                title="Oxford ST-expected IC vs Tree-assigned IC",
                row_levels=rl, col_levels=cl, counts=ct,
                max_rows=50, max_cols=50,
            )

    # Group col confusion matrices (cgMLST/BAPS/etc)
    if args.group_cols:
        cols = [c.strip() for c in args.group_cols.split(",") if c.strip()]
        for gc in cols:
            if gc not in fieldnames:
                logging.warning("Group column not found (skipping): %s", gc)
                continue

            mask_fn = lambda rr: (rr.get(gc, "") or "").strip() != "" and (rr.get(out_ic_col, "") or "").strip() != args.ua_label
            rl, cl, ct = crosstab(rows, gc, out_ic_col, mask_fn=mask_fn)

            safe_name = "".join([ch if ch.isalnum() or ch in ["_", "-"] else "_" for ch in gc])
            out_csv = os.path.join(args.outdir, f"confusion_{safe_name}_vs_assignedIC.csv")
            write_crosstab_csv(out_csv, rl, cl, ct, gc, "Assigned_IC (tree)")

            if args.plots:
                plot_confusion_heatmap(
                    out_png=os.path.join(args.outdir, f"confusion_{safe_name}_vs_assignedIC.png"),
                    out_pdf=os.path.join(args.outdir, f"confusion_{safe_name}_vs_assignedIC.pdf"),
                    title=f"{gc} vs Tree-assigned IC (counts)",
                    row_levels=rl, col_levels=cl, counts=ct,
                    max_rows=args.max_groups_plot, max_cols=min(args.max_groups_plot, 50),
                )

    # Write outputs
    safe_makedirs(args.outdir)
    write_summary(rows, args.ic_col, out_ic_col, args.ua_label, args.outdir)

    # Ensure output fieldnames include new cols (append if missing)
    out_fields = list(fieldnames)
    for c in extra_cols:
        if c not in out_fields:
            out_fields.append(c)

    logging.info("Writing output CSV: %s", args.output)
    write_csv_rows(args.output, rows, out_fields)

    # Also write an enriched copy into outdir for convenience
    enriched = os.path.join(args.outdir, os.path.basename(args.output).replace(".csv", "_enriched.csv"))
    write_csv_rows(enriched, rows, out_fields)
    logging.info("Wrote enriched output: %s", enriched)

    # Plots
    if args.plots:
        logging.info("Generating plots into: %s", args.outdir)
        plot_ic_before_after(rows, args.ic_col, out_ic_col, args.outdir)
        plot_histograms(rows, args.ic_col, out_ic_col, args.ua_label, args.outdir)

    logging.info("Done.")


if __name__ == "__main__":
    main()
