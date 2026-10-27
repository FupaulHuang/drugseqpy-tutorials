# Prepare and validate your inputs

[Tutorial home](index.md) · [Getting started](getting-started.md) · [QC](qc-workflow.md) · [Mixed species](mix-species-workflow.md) · [Single plate](single-plate-comparison.md) · [Multiple plates](multi-plate-comparison.md) · [API](api-reference.md)

Every workflow starts with raw UMI counts and a table that identifies each measured well. Complete these checks once, then continue to the [QC workflow](qc-workflow.md). The names **Sample 1–4** are neutral tutorial labels. The examples below describe a file contract; they do not include an experimental dataset.

## 1. Keep the matrix and metadata together

Use one directory per sample, with the three compressed matrix files directly inside its `matrix/` subdirectory:

```text
data/
├── sample_1/
│   ├── matrix/
│   │   ├── matrix.mtx.gz
│   │   ├── features.tsv.gz
│   │   └── barcodes.tsv.gz
│   └── metadata.csv
├── sample_2/
│   ├── matrix/
│   └── metadata.csv
├── sample_3/
│   ├── matrix/
│   └── metadata.csv
└── sample_4/
    ├── matrix/
    └── metadata.csv
```

The matrix bundle stores **features × wells**. Feature rows follow `features.tsv.gz`; well columns follow `barcodes.tsv.gz`. The loader returns a `DrugSeqData` object whose underlying AnnData uses **wells × features**, and keeps raw counts in `layers["counts"]`.

| Input | Contract |
| --- | --- |
| `matrix.mtx.gz` | Nonnegative, finite, integer UMI counts. Use raw counts, not CPM, normalized expression, log values, or z-scores. |
| `features.tsv.gz` | Stable feature ID in column 1; gene symbol in column 2 when available. Feature IDs should be unique. Symbols may repeat and are unsuitable as universal join keys. |
| `barcodes.tsv.gz` | One barcode per matrix column, in exactly the same order. |
| Metadata | One row per intended measured well. Join by an explicit barcode column; row order alone is insufficient. |

Record whether the matrix contains all assigned wells or only wells retained by an upstream filter. A well absent from a filtered matrix is not a measured zero. Keep the inclusion rule consistent when comparing samples and retain an exclusion table when subsetting.

## 2. Describe the experimental design explicitly

| Metadata field | Meaning and checks |
| --- | --- |
| `plate_id` | Physical plate identifier. Use `Sample 1`–`Sample 4` for the tutorial examples. Use distinct IDs for different plates; a sequencing sample and a physical plate are different concepts. |
| `barcode` | Matrix join key. It must be unique within the matrix being loaded. Preserve it as text. |
| `well_id` | Well label. If used for plotting, it must encode the correct physical position; see the coordinate check below. |
| `compound` | Exact condition label. Use consistent spelling and case for controls and treatments. |
| `dose`, `dose_unit`, `time` | Explicit experimental values when known. Keep unknown values missing; do not invent a dose or time from a label. |
| `sample_type` | Set explicitly where possible: `treatment`, `DMSO`, `vehicle`, `positive_ctrl`, `negative_ctrl`, or `media`. |
| `total_well` | Nominal plate capacity, such as 96; this is not the number of measured wells. |
| `x`, `y` | In this tutorial's input convention, `x` is physical column and `y` is physical row, both starting at 1. |
| `expected_species` | Optional designed species assignment for the mixed-species workflow. It should come from the experimental design. |

The loader can derive `sample_type` from a few exact `compound` labels, but that convenience does not establish the control design. Keep an explicit list of control well IDs, test conditions, exclusions, and the experimental unit that constitutes an independent replicate.

For multiple plates, include separate donor, batch, dose, and time fields where they affect the design. Repeated treatment labels within a plate are legitimate replicates; repeated barcode keys within one input matrix are ambiguous.

## 3. Check identity before importing

Run this from your analysis directory after installing `drugseqpy`. Every matrix well must have metadata. Metadata can also describe wells absent from the supplied matrix, for example blank wells omitted upstream; the example saves those rows in an inclusion audit. If exact equality is required by your experiment, enable the stricter check below. Keep unoccupied physical positions without barcodes in a separate layout annotation table.

```python
from pathlib import Path
import gzip

import numpy as np
import pandas as pd
from scipy import sparse
import drugseqpy as ds

matrix_dir = Path("data/sample_1/matrix")
metadata_path = Path("data/sample_1/metadata.csv")

metadata = pd.read_csv(metadata_path, dtype=str, keep_default_na=False)
required = {"plate_id", "barcode", "well_id", "compound", "sample_type"}
assert required <= set(metadata.columns), "Missing required tutorial metadata"
assert not metadata[list(required)].eq("").any().any(), "Empty identifiers or labels"

with gzip.open(matrix_dir / "barcodes.tsv.gz", "rt") as handle:
    matrix_ids = pd.Index([line.strip() for line in handle])
metadata_ids = pd.Index(metadata["barcode"])
assert metadata_ids.equals(metadata_ids.str.strip()), "Trim metadata barcode whitespace"

# The loader removes a terminal '-1'; check uniqueness after that operation.
matrix_keys = matrix_ids.str.replace(r"-1$", "", regex=True)
metadata_keys = metadata_ids.str.replace(r"-1$", "", regex=True)
assert matrix_keys.is_unique and metadata_keys.is_unique, "Duplicate normalized barcodes"
assert not matrix_keys.isin([""]).any(), "Empty matrix barcode"
assert set(matrix_keys) <= set(metadata_keys), "Matrix wells lack metadata"
require_exact_well_set = False
if require_exact_well_set:
    assert set(matrix_keys) == set(metadata_keys), "Matrix and metadata wells differ"

Path("results/input-checks").mkdir(parents=True, exist_ok=True)
metadata_audit = metadata.copy()
metadata_audit["in_matrix"] = metadata_keys.isin(matrix_keys)
metadata_audit.to_csv("results/input-checks/sample_1_metadata_audit.csv", index=False)
print("Metadata wells absent from matrix:", int((~metadata_audit["in_matrix"]).sum()))

with gzip.open(matrix_dir / "features.tsv.gz", "rt") as handle:
    features = pd.read_csv(handle, sep="\t", header=None, dtype=str)
assert features[0].notna().all() and features[0].is_unique, "Invalid feature IDs"

dsd = ds.create_drugseq_from_files(
    matrix_dir,
    metadata_path,
    metadata_key="barcode",
    min_overlap=1.0,
)
assert dsd.shape == (len(matrix_keys), len(features))
assert dsd.obs_names.equals(matrix_keys), "Unexpected observation order"
assert dsd.var_names.equals(pd.Index(features[0])), "Unexpected feature order"

counts = dsd.layers["counts"]
values = counts.data if sparse.issparse(counts) else np.asarray(counts).ravel()
assert np.isfinite(values).all(), "Nonfinite counts"
assert (values >= 0).all(), "Negative counts"
assert np.equal(values, np.floor(values)).all(), "Noninteger counts"

# Retain an independently calculated reference for subsequent QC tables.
raw_totals = np.asarray(counts.sum(axis=1)).ravel()
raw_detected = np.asarray((counts > 0).sum(axis=1)).ravel()
input_audit = dsd.obs[["plate_id", "barcode", "well_id", "compound"]].copy()
input_audit["raw_umi"] = raw_totals
input_audit["raw_detected_features"] = raw_detected
input_audit.to_csv("results/input-checks/sample_1_input_audit.csv")
```

`min_overlap=1.0` requires every matrix barcode to appear in metadata. It does not reject extra metadata rows. The explicit audit preserves those absent wells, while the optional equality assertion supports experiments that require exactly matching sets. The loader aligns rows by the join key rather than by the original CSV order.

These checks catch common format and alignment mistakes. Integer-valued expression alone cannot prove that a matrix contains raw counts; confirm how the matrix was produced as part of your private analysis record.

<a id="physical-coordinates"></a>

## 4. Verify physical plate positions

The current plate helper gives **parseable `well_id` values priority over `x/y`**. `A1` means row A, column 1. Numeric and `Well-dT` labels use column-major numbering: for an 8-row plate, 1–8 occupy the first column. A barcode label that looks like a numbered well can therefore produce a plausible but incorrect map when numbering differs from physical position.

Use the current layout annotations to reconcile positions. For an 8 × 12 plate, this example preserves the original well label and gives a plotting copy explicit physical coordinates:

```python
plot_dsd = dsd.copy()
plot_dsd.obs["well_label"] = plot_dsd.obs["well_id"].astype(str)

x = pd.to_numeric(plot_dsd.obs["x"], errors="raise")
y = pd.to_numeric(plot_dsd.obs["y"], errors="raise")
assert x.notna().all() and y.notna().all()
assert x.eq(x.round()).all() and y.eq(y.round()).all()
assert x.between(1, 12).all() and y.between(1, 8).all()
assert not pd.DataFrame({
    "plate_id": plot_dsd.obs["plate_id"], "x": x, "y": y
}).duplicated().any(), "Two observations occupy one physical position"

plot_dsd.obs["well_id"] = [
    f"{chr(64 + int(row))}{int(column)}" for row, column in zip(y, x)
]
layout = ds.build_plate_layout(plot_dsd.obs, nrow=8, ncol=12)
assert np.array_equal(layout["row"].to_numpy(), y.to_numpy())
assert np.array_equal(layout["column"].to_numpy(), x.to_numpy())
layout.to_csv("results/input-checks/sample_1_layout_audit.csv", index=False)
```

Use `plot_dsd` for plate figures and `label_col="well_label"` where the plotting function accepts it to display the original labels. Keep the analysis object and barcode identities unchanged. For another plate shape, adapt the row labels and bounds explicitly. Check a known occupied position and a known empty position against your design before interpreting spatial patterns. An empty position is different from an occupied well with zero counts.

## 5. Make observation IDs unique before merging

Different plates can reuse the same barcode sequence. Keep `barcode` as the original join key and give observations globally unique names before calling `merge_drugseq_objects`:

```python
# Run after loading each object; use sample_2, sample_3, or sample_4 as appropriate.
sample_key = "sample_1"
dsd.adata.obs_names = pd.Index([
    f"{sample_key}|{plate}|{barcode}"
    for plate, barcode in zip(dsd.obs["plate_id"], dsd.obs["barcode"])
])
dsd.obs["sample_id"] = dsd.obs_names
assert dsd.obs_names.is_unique
```

The merge function rejects duplicate observation IDs. It uses the union of features and fills missing features with zeros, so verify that the samples used compatible gene annotations before merging. A gene absent from a reference is not evidence that it was biologically unexpressed. For a controlled comparison, explicitly reconcile feature IDs and inspect the shared feature set first.

<a id="repeatability-and-concordance"></a>

## 6. Use terminology that matches the experimental design

| Term | Meaning and use in these tutorials |
| --- | --- |
| **Within-plate repeatability** | Descriptive consistency of matched wells measured under the same documented plate conditions. Profile correlations are one diagnostic, not a complete estimate of measurement precision. |
| **Between-plate consistency within one batch** | The scope of the measured Drug-seq case: matched conditions on three plates sharing a cell-preparation batch and repeated spatial layout. |
| **Cross-plate treatment-effect concordance** | Agreement of separately fitted treatment-versus-control effect vectors, assessed with correlation, direction and magnitude, and DEG overlap. It does not establish which experimental conditions were varied. |
| **Experimental reproducibility** | Consistency under specified changed measurement or experimental conditions. The claim must identify the conditions varied and the independent experimental units. |
| **Intermediate precision** | Within-laboratory precision across an extended period and specified changes such as runs, operators, or calibrations. Do not infer this solely from the presence of several plates. |
| **Computational reproducibility** | Obtaining consistent numerical results from the same inputs, analysis code, configuration, seeds, and compatible software environment. This does not establish experimental reproducibility. |

The measurement definitions distinguish repeatability under a specified stable setup, intermediate precision within a laboratory, and reproducibility under documented changed conditions. See the [VIM repeatability definition](https://jcgm.bipm.org/vim/en/2.20.html), [intermediate-precision definition](https://jcgm.bipm.org/vim/en/2.22.html), and [reproducibility definition](https://jcgm.bipm.org/vim/en/2.24.html).

For the measured case, the three plates share one cell-preparation batch, and biological independence of individual wells is not established. We therefore describe **within-plate repeatability** and **between-plate treatment-effect consistency within one batch**. Separately fitting models on those plates does not turn them into independent biological preparations. High correlation alone does not establish equal measurements, treatment causality, or reproducibility under changed conditions.

The package API names `compute_cross_plate_reproducibility()` and `plot_cross_plate_reproducibility()` are retained so code remains runnable. In tutorial prose their outputs are called *cross-plate treatment-effect concordance*. Reference signatures also retain the real argument names. Statements about reproducible code, seeds, or numerical results refer to computational reproducibility.

<a id="6-keep-comparisons-and-reports-reproducible"></a>

## 7. Keep analysis inputs and computational results reproducible

- **Record inclusion separately from QC flags.** Save the initial well set, upstream filtered status, your quality rule, and each exclusion reason. Write `AND` or `OR` explicitly when combining UMI and detected-gene cutoffs. A flag alone does not authorize exclusion.
- **State what each count measures.** Input read pairs, barcode-matched reads, deduplicated UMIs, and detected features are different quantities. Save the denominator of every percentage. If reporting only symbol-annotated genes, label that subset separately from all detected features.
- **Separate depth adjustment from normalization.** Equal-read comparisons require defined sampling and alignment inputs. Normalizing a count matrix does not recreate an equal-read sequencing experiment. For paired wells, record the matching key and report wells excluded from the intersection.
- **Keep expected species independent of the outcome.** A species inferred from a baseline sample is a derived label, not an independent designed truth. Record that distinction before interpreting agreement or contamination.
- **Preserve statistical decisions.** Save the exact control IDs, replicate unit, model formula, contrasts, filtering rule, normalization settings, gene-set snapshot, and tested-gene background. Pairwise correlations sharing a well are not independent replicates.
- **Save outputs that explain the figures.** Keep per-well QC and layout tables, complete contrast results, gene-resolution audits, software versions, and seeds with your analysis. Rebuilding HTML from saved results does not refit a model.

Before sharing a report, verify that its displayed totals agree with the saved tables, its figures load, its links resolve, and its controls and plate labels match the intended analysis. Test the navigation and tables at both desktop and narrow widths. A combined plate viewer changes presentation; it does not justify pooling otherwise incompatible designs.

Continue to [QC](qc-workflow.md), then choose [mixed-species QC](mix-species-workflow.md), [well-to-well repeatability](well-to-well-repeatability.md), [single-plate comparison](single-plate-comparison.md), or [multi-plate comparison](multi-plate-comparison.md). Use [bulk RNA-seq comparison](bulk-rnaseq-comparison.md) for independent-reference agreement. Use the [API reference](api-reference.md) to check parameters and return values.
