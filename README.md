# ICassigner.py

Conservative phylogeny-guided assignment of *Acinetobacter baumannii* **International Clones (ICs)** using a **core genome phylogeny** plus a metadata CSV containing (partial) IC labels.

Designed for harmonising IC labels for analysis/visualisation when:
- cgMLST is preferred for lineage definition, but
- papers frequently summarise lineages using **IC1/IC2/…**, and
- some isolates are labelled `UA` (unassigned).

## What it does

- **Never overwrites** existing IC labels (unless `--force_overwrite`)
- For `UA` isolates, walks up the tree to the **nearest ancestor clade** with sufficient labelled descendants.
- Assigns the **majority IC** only if conservative thresholds are met.

Default thresholds:
- `support_n >= 10` labelled descendants
- `majority_prop >= 0.90`

Otherwise isolates remain `UA`.

## Requirements

Minimal:
```bash
python -m pip install ete3
```

If you want plots (`--plots`):
```bash
python -m pip install matplotlib numpy
```

## Required inputs
### Core genome tree (Newick)

Tips must match sample IDs in metadata

Example:

```text
Acinetobacter-coreML.nwk
```

### Metadata CSV

Must contain:
- a sample ID column (matching tree tips)
- an IC column with known labels or `UA`

Example:

```csv
sample_id,IC,ST,Country
AB001,IC2,ST2,Russia
AB002,UA,ST78,Ukraine
```

### Basic usage (recommended)
```bash
python ICassigner.py \
  --tree Acinetobacter-coreML.nwk \
  --metadata metadata.csv \
  --tip_col sample_id \
  --ic_col IC \
  --plots \
  --outdir outputs_icassigner
```

Outputs:
- `metadata_with_conservative_IC.csv` (default)
- `outputs_icassigner/ICassigner_summary.txt`
- plots (PNG+PDF) if `--plots`

### Adjusting stringency (optional)

More conservative (publication quality):

```bash
--min_support 15 --min_prop 0.95
```

Slightly more permissive (exploratory analyses):

```bash
--min_support 5 --min_prop 0.8
```

### Optional: MLST cross-checks (Pasteur/Oxford)

If your metadata includes Pasteur and/or Oxford STs, ICassigner can:
- compute expected IC from canonical ST→IC anchors (sanity check)
- write conflict flags
- write confusion matrix CSVs (+ heatmap plots if `--plots`)

Example:

```bash
python ICassigner.py \
  --tree RAxML-result.Acinetobacter-coreML.nwk \
  --metadata FullMicroreactWyr-with-Russian-Metadata.csv \
  --tip_col sample_id --ic_col IC \
  --pasteur_st_col ST_Pasteur \
  --oxford_st_col ST_Oxford \
  --plots --outdir outputs_icassigner
```

### Optional: group cross-checks (cgMLST / BAPS / other)

Compare any grouping variable(s) to assigned ICs using `--group_cols`.

Example:

```bash
python ICassigner.py \
  --tree RAxML-result.Acinetobacter-coreML.nwk \
  --metadata FullMicroreactWyr-with-Russian-Metadata.csv \
  --tip_col sample_id --ic_col IC \
  --group_cols cgMLST_group,BAPS \
  --plots --outdir outputs_icassigner
```

Large group columns can produce huge matrices; ICassigner always writes full CSVs,
but plots are truncated to the top categories (default 30; configurable via `--max_groups_plot`).

### Output columns added

| Column                       | Meaning                                                 |
| ---------------------------- | ------------------------------------------------------- |
| `IC_tree_conservative`       | Final IC call (original IC or inferred IC)              |
| `IC_tree_support_n`          | Number of labelled descendants used                     |
| `IC_tree_majority_n`         | Count supporting the majority IC                        |
| `IC_tree_support_prop`       | Majority proportion supporting assigned IC              |
| `IC_tree_node_size`          | Size of the ancestor clade used (number of tips)        |
| `IC_expected_from_PasteurST` | Expected IC based on Pasteur ST anchors (if provided)   |
| `IC_PasteurST_conflict`      | `1` if expected != assigned (assigned not UA), else `0` |
| `IC_expected_from_OxfordST`  | Expected IC based on Oxford ST anchors (if provided)    |
| `IC_OxfordST_conflict`       | `1` if expected != assigned (assigned not UA), else `0` |

If an isolate remains `UA`, this is by design and indicates genuine phylogenetic ambiguity.

### Outputs written to `--outdir`
- `ICassigner_summary.txt`
- `tips_missing_in_tree.txt` (if applicable)
- confusion matrix CSVs (e.g., `confusion_expectedIC_PasteurST_vs_assignedIC.csv`)
- plots (if `--plots`):
  - `IC_counts_before_after.png/.pdf`
  - `IC_inference_support_hist.png/.pdf`
  - `IC_inference_prop_hist.png/.pdf`
  - `confusion_*.png/.pdf`

## Suggested citation text (software)

International clone assignment was performed using **ICassigner.py**, a conservative phylogeny-guided IC propagation tool based on majority support within the core genome tree.  

First described in *Pascoe & Mourkas et al.* (in preparation).
