# Getting started

[Learning path](index.md) · [Next: prepare inputs](input-preparation.md)

Use this chapter to get a working environment and run the synthetic practice workflow. You need a terminal, basic Python familiarity, and Python 3.10 or later. The examples were checked with drugseqpy 0.1.0. Feature availability matters because different development revisions can report the same version.

## 1. Install drugseqpy

Download or clone this tutorial repository, open a terminal in its top-level folder, then create an environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install "git+https://github.com/FupaulHuang/drugseqpy.git"
```

On Windows, activate with `.venv\Scripts\activate`. The install URL is the [software repository](https://github.com/FupaulHuang/drugseqpy). Python packages used by the examples are installed with drugseqpy. Excel input additionally needs an appropriate pandas engine, such as `openpyxl` for `.xlsx`; the bundled metadata are CSV.

Check the functions required by this tutorial before starting:

```bash
python scripts/check_environment.py
```

These tutorials were validated against a local development build; a matching public release is not assumed. If a function is missing, update to a package revision that exports it, or install a compatible checkout with `python -m pip install /path/to/drugseqpy`. A version string of `0.1.0` alone is insufficient to establish compatibility. Save `python -m pip freeze` with your results to record the exact environment. The default examples use the package's Python `limma_voom` implementation; the alternative edgeR differential-expression backend needs R, edgeR, and rpy2.

## 2. Inspect the practice inputs

Four small, entirely synthetic matrix bundles are included under `data/`. Each contains 48 wells and 800 features. A matrix is stored as genes × wells; drugseqpy imports it into an object whose rows are wells and columns are genes. Each `metadata.csv` joins to the corresponding matrix through `barcode`.

```text
data/
├── sample_1/metadata.csv
├── sample_1/matrix/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}
├── sample_2/…
├── sample_3/…
├── sample_4/…
└── practice-manifest.json
```

Sample 1 and Sample 2 have matching DMSO, Drug1, Drug2, and water groups. Sample 3 and Sample 4 have `231` (human), `3T3` (mouse), and water labels. These labels demonstrate how metadata drive the analysis. Gene IDs and counts are synthetic; marker symbols provide recognizable plotting examples. Do not use these data to draw biological conclusions.

To generate another practice copy without overwriting your inputs:

```bash
python scripts/make_practice_data.py --output practice-copy
```

## 3. Run QC, species review, and comparison

Run from the tutorial repository root. Start with QC:

```bash
python scripts/run_tutorial.py qc --data-dir data --output-dir results
```

Then run the other two workflows:

```bash
python scripts/run_tutorial.py species --data-dir data --output-dir results
python scripts/run_tutorial.py comparison --data-dir data --output-dir results
```

Or run all three with `python scripts/run_tutorial.py all`. The runner writes images and tables to `results/`. It does not change the input matrices. Open the PNGs in an image viewer and the CSV files in a spreadsheet or Python.

| Workflow | Key outputs | First review |
| --- | --- | --- |
| QC | `results/qc/sample_qc.csv`, QC distributions and plate maps | Which wells have low depth, few detected genes, or unusual positions? |
| Species | `results/species/species_metrics.csv`, species figures | Do intended and observed species agree, and are fractions supported by enough counts? |
| Comparison | `results/comparison/design_counts.csv`, `batch_metrics.csv`, DE outputs | Are treatment and control replicated on both plates, and are effects consistent? |

The runner is a concise starting workflow. The detailed chapters add more figure variants and explain the archived worked results. QC uses Sample 3 and Sample 4; comparison computes its own QC on Sample 1 and Sample 2. Existing figure/result names may be overwritten when you rerun into the same output directory, so use a new output directory for each analysis you want to keep.

## 4. Read the results before changing settings

The practice runner uses a dominant-species threshold of 0.9 and a minimum of 100 species-assigned UMIs. These are illustrative review settings. The archived mixed-species figures use 0.5 and 0 UMIs, so their call categories differ. Fractions and raw counts remain useful at either setting; call labels must be interpreted with their thresholds.

The runner flags unusual DMSO controls without automatically excluding them. It excludes blank groups from treatment models and requires at least two wells per treatment and DMSO on each plate. A technical replicate is not automatically an independent biological replicate. Confirm the experimental unit before inference.

For reviewed exclusions, optional embeddings, local GMT enrichment, species mappings, and plate dimensions, inspect the explicit options:

```bash
python scripts/run_tutorial.py --help
```

The default comparison includes PCA, repeatability, cross-plate reproducibility, and differential expression. Add `--full-embeddings` for UMAP and t-SNE. Enrichment and STARsolo read-level statistics are not required for the first successful run.

## 5. Bring your own inputs

Use a separate data directory with the same folder structure. Review [Prepare inputs](input-preparation.md) first, particularly barcode matching and physical coordinates. The runner relabels the four datasets as Sample 1–4; it is tailored to these example roles and labels. The chapter code and API reference show how to generalize to other controls, species, plate counts, doses, and time points.

[← Learning path](index.md) · [Prepare inputs →](input-preparation.md)
