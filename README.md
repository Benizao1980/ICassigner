# How to use ICassigner.py

## Requirements

```bash
pip install pandas ete3
```

(ETE3 does not require graphical dependencies for tree parsing.)

## Required inputs
### Core genome tree (Newick)

Tips must match sample IDs in metadata

Example:

```text
RAxML-result.Acinetobacter-coreML.nwk
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
  --tree RAxML-result.Acinetobacter-coreML.nwk \
  --metadata FullMicroreactWyr-with-Russian-Metadata.csv
```

This will:
- keep all existing IC labels
- conservatively assign ICs to UA isolates
- write:
```text
metadata_with_conservative_IC.csv
```

### Adjusting stringency (optional)

More conservative (publication quality):

```bash
--min_support 15 --min_prop 0.95
```

Slightly more permissive (exploratory analyses):

```
--min_support 5 --min_prop 0.8
```

### Output columns added

| Column                  | Meaning                                               |
|-------------------------|-------------------------------------------------------|
| IC_tree_conservative    | Final IC call (original IC or inferred IC)            |
| IC_tree_support_n       | Number of labelled neighbours used                    |
| IC_tree_support_prop    | Proportion supporting assigned IC                    |


If an isolate remains `UA`, this is by design and indicates genuine phylogenetic ambiguity.

## Suggested citation text (software)

International clone assignment was performed using **ICassigner.py**, a conservative phylogeny-guided IC propagation tool based on majority support within the core genome tree.  

First described in *Pascoe & Mourkas et al.* (in preparation).
