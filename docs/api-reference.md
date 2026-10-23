# API reference

Use this reference alongside the [learning path](index.md), [QC](qc-workflow.md), [mixed-species](mix-species-workflow.md), [single-plate](single-plate-comparison.md), and [multi-plate](multi-plate-comparison.md) chapters. Functions are grouped by task. Each entry documents the reviewed Python signature, default values, parameters, outputs, side effects, and workflow-specific notes. Signature blocks describe the API; use the workflow chapters for executable calls.

**Compatibility.** Run `python scripts/check_environment.py` before starting. These contracts were checked against a development build reporting version 0.1.0; a matching public release is not assumed. Signatures below are reference displays, not complete runnable cells.

**Import convention.** Use `import drugseqpy as ds`, then call any public function directly—for example, `ds.merge_drugseq_objects()` or `ds.generate_qc_report()`. The names below are all exported from the package root.

## Read the data object

`DrugSeqData` wraps an AnnData object. `dsd.adata.obs` contains wells, `dsd.adata.var` contains genes, and `dsd.adata.layers['counts']` retains raw UMI counts. `dsd.adata.X` starts as counts on import and changes when you normalize or correct a view. `dsd.copy()` creates an independent copy; `dsd.subset(obs_mask=..., var_mask=...)` selects wells or genes. Set observation names through `dsd.adata.obs_names` when making them globally unique before merging.

Follow [Getting started](getting-started.md) for the executable workflow. The reference signatures below describe arguments and defaults; use the chapter examples for complete runnable cells.

<a id="input-objects"></a>

## Input and data objects

Read expression and metadata, construct DrugSeqData, and merge plates.

| Function | Description |
| --- | --- |
| [`drugseqpy.load_expression_matrix()`](#api-load-expression-matrix) | Read either a gene-by-sample table or a STARsolo count-matrix bundle. |
| [`drugseqpy.load_metadata()`](#api-load-metadata) | Read sample metadata from CSV, TSV/TXT, or Excel. |
| [`drugseqpy.create_drugseq_from_files()`](#api-create-drugseq-from-files) | Load, normalize barcode identifiers, align counts to metadata, and construct a DrugSeqData object. |
| [`drugseqpy.activate_gene_names()`](#api-activate-gene-names) | Activate a gene-annotation column as the AnnData variable names without changing matrix values or gene order. |
| [`drugseqpy.merge_drugseq_objects()`](#api-merge-drugseq-objects) | Combine separately loaded plates along the sample axis while retaining the union of genes. |

<a id="input-objects-details"></a>

### Detailed reference

<a id="api-load-expression-matrix"></a>

### `drugseqpy.load_expression_matrix()`

Read either a gene-by-sample table or a STARsolo count-matrix bundle.

```python
load_expression_matrix(path: 'str | Path') -> 'tuple[pd.DataFrame, pd.DataFrame | None]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `path` | str \| Path | required | File path to a CSV/TSV count table, or a directory containing `matrix.mtx.gz`, `features.tsv.gz`, and `barcodes.tsv.gz`. |

**Returns / side effects:** A gene-by-sample sparse DataFrame and optional gene annotation table. Barcode suffix `-1` is removed.

Point this loader at the directory containing the three matrix files. Both raw and filtered STARsolo bundles use this format; record which well set you selected. The whole STAR run root is used separately for read-level QC.

<a id="api-load-metadata"></a>

### `drugseqpy.load_metadata()`

Read sample metadata from CSV, TSV/TXT, or Excel.

```python
load_metadata(path: 'str | Path', index_col: 'str | int | None' = None, delimiter: 'str | None' = 'auto') -> 'pd.DataFrame'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `path` | str \| Path | required | Path to `.csv`, `.tsv`, `.txt`, `.xlsx`, or `.xls` metadata. |
| `index_col` | str \| int \| None | None | Column name or zero-based column number used as the sample index. If omitted, `sample_id` is preferred, then `barcode`. |
| `delimiter` | str \| None | 'auto' | Metadata delimiter. `'auto'` uses comma for CSV, tab for TSV, and delimiter detection for TXT. An explicit value overrides the extension for text files; Excel does not accept this parameter. |

**Returns / side effects:** A pandas DataFrame indexed by the selected sample identifier; duplicate identifiers raise an error.

<a id="api-create-drugseq-from-files"></a>

### `drugseqpy.create_drugseq_from_files()`

Load, normalize barcode identifiers, align counts to metadata, and construct a DrugSeqData object.

```python
create_drugseq_from_files(expression_path: 'str | Path', metadata_path: 'str | Path', *, metadata_key: 'str | int | None' = None, delimiter: 'str | None' = 'auto', min_overlap: 'float' = 0.5) -> 'DrugSeqData'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `expression_path` | str \| Path | required | Count table or a directory containing `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`. |
| `metadata_path` | str \| Path | required | CSV, TSV/TXT, or Excel metadata path. |
| `metadata_key` | str \| int \| None | None | Metadata column used for count/metadata matching. `None` uses `sample_id` or `barcode` automatically. |
| `delimiter` | str \| None | 'auto' | Delimiter used for the metadata file. `'auto'` detects TXT delimiters; an explicit separator overrides extension-based parsing for delimited text metadata. |
| `min_overlap` | float | 0.5 | Minimum fraction of expression barcodes that must match metadata. A lower overlap raises `ValueError`. |

**Returns / side effects:** A DrugSeqData object with samples in `adata.obs`, genes in `adata.var`, and raw counts in `adata.layers['counts']`. It also creates `sample_id = plate_id + barcode` and `treatment = compound + dose + time`, skipping empty components.

<a id="api-activate-gene-names"></a>

### `drugseqpy.activate_gene_names()`

Activate a gene-annotation column as the AnnData variable names without changing matrix values or gene order.

```python
activate_gene_names(dsd: 'DrugSeqData', label_col: 'str' = 'gene_symbol', *, stable_id_col: 'str' = 'ensembl_id', inplace: 'bool' = True) -> 'DrugSeqData | None'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData whose active gene labels will be changed. |
| `label_col` | str | 'gene\_symbol' | Column in `dsd.var` to activate, such as `gene_symbol`. Use the stable-ID column to switch back. |
| `stable_id_col` | str | 'ensembl\_id' | Column that preserves the original variable names and supplies fallbacks for missing labels. |
| `inplace` | bool | True | Modify the input and return `None`; `False` returns an independently modified copy. |

**Returns / side effects:** A modified DrugSeqData copy when `inplace=False`; otherwise `None`. Duplicate active labels receive unique suffixes and the active-label configuration is recorded in `uns['active_gene_names']`.

<a id="api-merge-drugseq-objects"></a>

### `drugseqpy.merge_drugseq_objects()`

Combine separately loaded plates along the sample axis while retaining the union of genes.

```python
merge_drugseq_objects(*objects: 'DrugSeqData', fill_missing: 'float' = 0.0, batch_key: 'str' = 'batch') -> 'DrugSeqData'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `objects` | DrugSeqData | one or more | Two or more DrugSeqData objects. Sample names must already be unique across plates. |
| `fill_missing` | float | 0.0 | Fill value for features absent from an input object. Confirm annotation compatibility: an absent feature is not evidence of zero biological expression. |
| `batch_key` | str | 'batch' | Name of the new `obs` column recording which input object supplied each sample. |

**Returns / side effects:** A new merged DrugSeqData object. Shared gene annotations are preserved when feature indexes match.

<a id="quality-control"></a>

## Quality-control metrics

Calculate and summarize sample-, metadata-, and plate-level QC information.

| Function | Description |
| --- | --- |
| [`drugseqpy.compute_qc_metrics()`](#api-compute-qc-metrics) | Calculate sample-level library composition, housekeeping stability, Gini inequality, and a multivariate MAD outlier score. |
| [`drugseqpy.validate_metadata()`](#api-validate-metadata) | Check required columns, missing values, and recognized sample-type labels. |
| [`drugseqpy.compute_plate_qc()`](#api-compute-plate-qc) | Summarize control separation and variability within each plate. |
| [`drugseqpy.plate_qc_summary()`](#api-plate-qc-summary) | Convert stored plate-QC results to a tidy table. |

<a id="quality-control-details"></a>

### Detailed reference

<a id="api-compute-qc-metrics"></a>

### `drugseqpy.compute_qc_metrics()`

Calculate sample-level library composition, housekeeping stability, Gini inequality, and a multivariate MAD outlier score.

```python
compute_qc_metrics(dsd: 'DrugSeqData', mito_pattern: 'str | None' = None, ribo_pattern: 'str | None' = None, hk_genes: 'list[str] | None' = None, min_hk_detected: 'int' = 5, inplace: 'bool' = True, *, mito_prefixes: 'str | Sequence[str] | None' = ('MT-', 'mt-'), ribo_prefixes: 'str | Sequence[str] | None' = ('RPL', 'RPS', 'rpl', 'rps'), mito_genes: 'Sequence[str] | None' = None, ribo_genes: 'Sequence[str] | None' = None) -> 'DrugSeqData | None'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData containing raw counts in `layers['counts']` (preferred) or `X`. |
| `mito_pattern` | str \| None | None | Optional regular expression matched against both gene symbols and feature IDs. It overrides `mito_prefixes` unless `mito_genes` is supplied. |
| `ribo_pattern` | str \| None | None | Optional regular expression matched against both gene symbols and feature IDs. It overrides `ribo_prefixes` unless `ribo_genes` is supplied. |
| `hk_genes` | list[str] \| None | None | Optional housekeeping-gene symbols. `None` uses the package list/fallback list. |
| `min_hk_detected` | int | 5 | Minimum number of supplied housekeeping genes present before housekeeping metrics are calculated. |
| `inplace` | bool | True | If `True`, write metrics into the existing object; if `False`, calculate them on and return a copy. |
| `mito_prefixes` | str \| Sequence[str] \| None | ('MT-', 'mt-') | String or sequence of mitochondrial prefixes. The default recognizes both `MT-` and `mt-`. |
| `ribo_prefixes` | str \| Sequence[str] \| None | ('RPL', 'RPS', 'rpl', 'rps') | String or sequence of ribosomal prefixes. The default recognizes uppercase and lowercase `RPL`/`RPS`. |
| `mito_genes` | Sequence[str] \| None | None | Optional exact mitochondrial gene symbols or feature IDs. A supplied list has highest priority and disables prefix/regex selection for this set. |
| `ribo_genes` | Sequence[str] \| None | None | Optional exact ribosomal gene symbols or feature IDs. A supplied list has highest priority and disables prefix/regex selection for this set. |

**Returns / side effects:** Adds `total_umi`, `n_genes_det`, `pct_mito`, `pct_ribo`, housekeeping metrics, `gini_index`, and `outlier_score` to `obs`; marks `var['is_mito']`/`var['is_ribo']`; and records the resolved gene definitions in `uns['qc_gene_definitions']`. Returns only when `inplace=False`.

Selection priority is exact gene list, then regex, then prefixes. The outlier score is the Euclidean norm of five median-centered, unscaled-MAD-normalized QC metrics; no 1.4826 consistency multiplier is used.

<a id="api-validate-metadata"></a>

### `drugseqpy.validate_metadata()`

Check required columns, missing values, and recognized sample-type labels.

```python
validate_metadata(dsd: 'DrugSeqData', required_cols: 'list[str] | None' = None, verbose: 'bool' = True) -> 'dict'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData whose `obs` metadata will be checked. |
| `required_cols` | list[str] \| None | None | Columns that must exist and contain no missing values. `None` uses the package defaults. |
| `verbose` | bool | True | Print issues and per-plate sample-type counts when enabled. |

**Returns / side effects:** A dictionary with Boolean `ok` and a list of human-readable `issues`.

<a id="api-compute-plate-qc"></a>

### `drugseqpy.compute_plate_qc()`

Summarize control separation and variability within each plate.

```python
compute_plate_qc(dsd: 'DrugSeqData', signal_col: 'str' = 'log10_total_umi', pos_ctrl_label: 'str' = 'positive_ctrl', neg_ctrl_label: 'str' = 'DMSO', inplace: 'bool' = True) -> 'DrugSeqData | None'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData after `compute_qc_metrics()`. |
| `signal_col` | str | 'log10\_total\_umi' | Numeric `obs` column used as the assay signal; the default is log10 UMI. |
| `pos_ctrl_label` | str | 'positive\_ctrl' | Value in `sample_type` identifying positive-control wells. |
| `neg_ctrl_label` | str | 'DMSO' | Value in `sample_type` identifying negative/reference-control wells. |
| `inplace` | bool | True | Store results in the current object when `True`, otherwise return a modified copy. |

**Returns / side effects:** Stores per-plate Z-prime, SSMD, negative-control CV/correlation, spatial matrices, and pass/warn/fail flags in `adata.uns['plate_qc']`.

Z-prime and SSMD require suitable positive and negative controls. Sample 3/Sample 4 are species-mixture plates, so missing control classes produce `NaN` and a warning flag; this is expected and must not be treated as an assay failure.

<a id="api-plate-qc-summary"></a>

### `drugseqpy.plate_qc_summary()`

Convert stored plate-QC results to a tidy table.

```python
plate_qc_summary(dsd: 'DrugSeqData') -> 'pd.DataFrame'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData after `compute_plate_qc()`. |

**Returns / side effects:** A DataFrame with one row per plate; raises an error when plate QC has not been computed.

<a id="subsetting-filtering"></a>

## Subsetting and filtering

Optionally remove reviewed low-quality wells or low-expression genes.

| Function | Description |
| --- | --- |
| [`drugseqpy.filter_samples()`](#api-filter-samples) | Return a new object after applying reviewed sample/well QC thresholds. |
| [`drugseqpy.filter_genes()`](#api-filter-genes) | Return a new object after removing low-expression and optionally mitochondrial/ribosomal genes. |

<a id="subsetting-filtering-details"></a>

### Detailed reference

<a id="api-filter-samples"></a>

### `drugseqpy.filter_samples()`

Return a new object after applying reviewed sample/well QC thresholds.

```python
filter_samples(dsd: 'DrugSeqData', min_umi: 'float' = 5000, max_umi: 'float' = inf, min_genes: 'int' = 200, max_pct_mito: 'float' = 25.0, max_pct_ribo: 'float' = inf, min_hk_mean: 'float' = 0.0, max_hk_cv: 'float' = inf, max_outlier_score: 'float' = 10.0, keep_sample_types: 'list[str] | None' = None, verbose: 'bool' = True) -> 'DrugSeqData'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData with sample QC metrics already calculated by `compute_qc_metrics()`. |
| `min_umi` | float | 5000 | Minimum total UMI count retained per unprotected well. |
| `max_umi` | float | inf | Maximum total UMI count retained per unprotected well. |
| `min_genes` | int | 200 | Minimum detected-gene count retained per unprotected well. |
| `max_pct_mito` | float | 25.0 | Maximum mitochondrial UMI percentage. |
| `max_pct_ribo` | float | inf | Maximum ribosomal UMI percentage. |
| `min_hk_mean` | float | 0.0 | Minimum housekeeping-gene mean log expression when that metric is available. |
| `max_hk_cv` | float | inf | Maximum housekeeping-gene coefficient of variation when available. |
| `max_outlier_score` | float | 10.0 | Maximum multivariate MAD outlier score. |
| `keep_sample_types` | list[str] \| None | None | Sample types exempt from thresholds. `None` protects DMSO, vehicle, and positive/negative controls; pass an empty list only after deciding no wells require protection. |
| `verbose` | bool | True | Print retained/removed counts and failure reasons. |

**Returns / side effects:** A new filtered DrugSeqData object; the input object is unchanged.

Thresholds are experiment-specific. Preview the failing-well mask and protect required controls before applying this optional operation.

<a id="api-filter-genes"></a>

### `drugseqpy.filter_genes()`

Return a new object after removing low-expression and optionally mitochondrial/ribosomal genes.

```python
filter_genes(dsd: 'DrugSeqData', min_count: 'int' = 1, min_samples: 'int' = 2, group_aware: 'bool' = True, group_col: 'str' = 'compound', min_cpm: 'float' = 0.0, min_cpm_samples: 'int' = 2, remove_mito: 'bool' = False, remove_ribo: 'bool' = False, mito_pattern: 'str' = '^MT-', ribo_pattern: 'str' = '^RP[SL]', verbose: 'bool' = True) -> 'DrugSeqData'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData containing raw counts in `layers['counts']`. |
| `min_count` | int | 1 | Minimum raw count used to call a gene expressed in one sample. |
| `min_samples` | int | 2 | Minimum samples meeting `min_count` globally or within at least one group. |
| `group_aware` | bool | True | When `True`, retain genes passing in at least one group so treatment-induced genes are not lost. |
| `group_col` | str | 'compound' | Metadata column defining groups for group-aware filtering; defaults to `compound`. |
| `min_cpm` | float | 0.0 | Optional CPM threshold; zero disables CPM filtering. |
| `min_cpm_samples` | int | 2 | Minimum samples that must meet `min_cpm`. |
| `remove_mito` | bool | False | Optionally remove genes whose feature names match `mito_pattern`. |
| `remove_ribo` | bool | False | Optionally remove genes whose feature names match `ribo_pattern`. |
| `mito_pattern` | str | '^MT-' | Regular expression used against feature names when `remove_mito=True`. |
| `ribo_pattern` | str | '^RP[SL]' | Regular expression used against feature names when `remove_ribo=True`. |
| `verbose` | bool | True | Print retained and removed gene counts. |

**Returns / side effects:** A new gene-filtered DrugSeqData object; the input object is unchanged.

Group-aware filtering is usually safer for perturbation screens. Pattern-based mitochondrial/ribosomal removal requires feature names compatible with the supplied patterns.

<a id="plate-layouts"></a>

## Plate layouts

Resolve physical coordinates and generate spatial plate visualizations.

| Function | Description |
| --- | --- |
| [`drugseqpy.build_plate_layout()`](#api-build-plate-layout) | Resolve each metadata row to a physical plate coordinate and attach the requested plotting value. |
| [`drugseqpy.plot_plate_layout()`](#api-plot-plate-layout) | Draw one spatial panel per plate with adaptive marker, font, panel, legend, and colorbar sizing. |
| [`drugseqpy.generate_plate_layout_figures()`](#api-generate-plate-layout-figures) | Generate one universal plate-layout figure and save its resolved plotting table. |

<a id="plate-layouts-details"></a>

### Detailed reference

<a id="api-build-plate-layout"></a>

### `drugseqpy.build_plate_layout()`

Resolve each metadata row to a physical plate coordinate and attach the requested plotting value.

```python
build_plate_layout(obs: 'pd.DataFrame', value_col: 'str | None' = None, *, nrow: 'int' = 8, ncol: 'int' = 12) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `obs` | pd.DataFrame | required | Sample metadata DataFrame. Parseable `well_id` values take precedence over `x/y`; validate the resulting positions against the physical design. |
| `value_col` | str \| None | None | Optional metadata/QC column copied into the generic `value` field. |
| `nrow` | int | 8 | Number of physical plate rows. The default is 8 for a 96-well plate; use 16 for a standard 384-well plate. |
| `ncol` | int | 12 | Number of physical plate columns. The default is 12 for a 96-well plate; use 24 for a standard 384-well plate. |

**Returns / side effects:** A tidy layout DataFrame with numeric `row`/`column`, display labels, and missing-well status.

For numbered Well-dT identifiers, numbering is column-major: Well-dT-1 through Well-dT-8 run down column 1, then Well-dT-9 starts column 2.

<a id="api-plot-plate-layout"></a>

### `drugseqpy.plot_plate_layout()`

Draw one spatial panel per plate with adaptive marker, font, panel, legend, and colorbar sizing.

```python
plot_plate_layout(obs: 'pd.DataFrame', value_col: 'str | None' = None, *, well_color: 'str | None' = None, plate_id: 'str' = 'plate_id', split_by: 'str | None' = None, label_col: 'str | None' = None, continuous: 'bool' = False, log_values: 'bool' = True, show_legend: 'bool' = True, nrow: 'int' = 8, ncol: 'int' = 12, nrow_plate: 'int | None' = None, figsize: 'tuple[float, float] | None' = None, well_size: 'float' = 0.76, well_border: 'bool' = True, well_border_color: 'str' = '#111827', well_border_width: 'float' = 0.8, palette: 'str | list[str] | tuple[str, ...] | dict[str, str] | None' = None, cmap: 'str' = 'RdYlBu_r', label_color: 'str' = '#111827', label_size: 'float | None' = None, label_style: 'str' = 'normal', label_weight: 'str | int' = 'normal')
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `obs` | pd.DataFrame | required | Metadata plus the numeric or categorical field selected by `value_col`. |
| `value_col` | str \| None | None | Metadata/QC column controlling well color. Supply exactly one of `value_col` or `well_color`. |
| `well_color` | str \| None | None | Literal fixed fill color, such as `'white'`. This is not a metadata column and does not create a legend. |
| `plate_id` | str | 'plate\_id' | Name of the metadata column containing the real plate identities. Defaults to `'plate_id'`. |
| `split_by` | str \| None | None | Optional metadata column creating separate physical-layout panels for each plate-by-split combination. `None` uses one panel per plate. |
| `label_col` | str \| None | None | Optional column printed inside occupied wells. `None` is the default and draws no labels. |
| `continuous` | bool | False | Use a continuous colormap and colorbar when `True`. The default `False` treats `value_col` as categorical. |
| `log_values` | bool | True | Apply `log1p` to continuous values before coloring. Defaults to `True`; the source table retains the original value and the transformed plotting value. |
| `show_legend` | bool | True | Show the categorical `value_col` legend or the continuous colorbar. Fixed-color plots have neither. |
| `nrow` | int | 8 | Number of physical plate rows; defaults to 8. |
| `ncol` | int | 12 | Number of physical plate columns; defaults to 12. |
| `nrow_plate` | int \| None | None | Maximum number of plate panels in one figure row. For two plates, `1` stacks them vertically and `2` places them side by side; `None` uses the adaptive arrangement. |
| `figsize` | tuple[float, float] \| None | None | Optional explicit figure size in inches. `None` adapts to plate dimensions and panel count. |
| `well_size` | float | 0.76 | Well diameter relative to one grid cell. The adaptive default is 0.76; values up to 1.25 are accepted. |
| `well_border` | bool | True | Draw borders around every occupied and missing well. Defaults to `True`. |
| `well_border_color` | str | '#111827' | Matplotlib color used for well borders when enabled. |
| `well_border_width` | float | 0.8 | Well-border line width in points. |
| `palette` | str \| list[str] \| tuple[str, ...] \| dict[str, str] \| None | None | Categorical colors as a Seaborn palette name, a color sequence, or a mapping from category labels to colors. |
| `cmap` | str | 'RdYlBu\_r' | Matplotlib colormap used when `continuous=True`. The default `RdYlBu_r` runs from bright blue through yellow to red. |
| `label_color` | str | '#111827' | Matplotlib color for text drawn inside wells. |
| `label_size` | float \| None | None | Optional label font size in points. `None` adapts to well size and plate format. |
| `label_style` | str | 'normal' | Font style for labels, for example `'normal'` or `'italic'`. |
| `label_weight` | str \| int | 'normal' | Font weight for labels, for example `'normal'`, `'bold'`, or a numeric weight. |

**Returns / side effects:** A `(matplotlib.figure.Figure, layout_dataframe)` tuple. Missing wells are white and use the configured border settings.

<a id="api-generate-plate-layout-figures"></a>

### `drugseqpy.generate_plate_layout_figures()`

Generate one universal plate-layout figure and save its resolved plotting table.

```python
generate_plate_layout_figures(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), plate_value: 'str | None' = None, plate_gene: 'str | Sequence[str] | None' = None, gene_aggregation: 'str' = 'mean', well_color: 'str | None' = None, plate_id: 'str' = 'plate_id', split_by: 'str | None' = None, label_col: 'str | None' = None, continuous: 'bool' = False, log_values: 'bool' = True, show_legend: 'bool' = True, nrow: 'int' = 8, ncol: 'int' = 12, nrow_plate: 'int | None' = None, output_name: 'str | None' = None, figsize: 'tuple[float, float] | None' = None, well_size: 'float' = 0.76, well_border: 'bool' = True, well_border_color: 'str' = '#111827', well_border_width: 'float' = 0.8, palette: 'PaletteSpec' = None, cmap: 'str' = 'RdYlBu_r', label_color: 'str' = '#111827', label_size: 'float | None' = None, label_style: 'str' = 'normal', label_weight: 'str | int' = 'normal') -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData containing aligned metadata and raw counts. Missing QC metrics are computed automatically. |
| `output_dir` | str \| Path | required | Directory receiving the plate-layout image files and coordinate source tables. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `plate_value` | str \| None | None | Metadata/QC field controlling well color. Exactly one of `plate_value`, `plate_gene`, or `well_color` is required. |
| `plate_gene` | str \| Sequence[str] \| None | None | One gene symbol/feature ID, or a sequence defining a gene set. A string plots that gene alone; a sequence is aggregated per well. |
| `gene_aggregation` | str | 'mean' | Gene-set aggregation: `'mean'` (default), `'sum'`, or `'median'`. It is ignored for a single string gene. |
| `well_color` | str \| None | None | Literal fixed fill color, such as `'white'`; it is not a metadata column. |
| `plate_id` | str | 'plate\_id' | Name of the metadata column containing actual plate identities. Defaults to `'plate_id'`. |
| `split_by` | str \| None | None | Optional metadata column creating separate plate-by-split layout panels. `None` keeps one panel per physical plate. |
| `label_col` | str \| None | None | Optional metadata field printed in occupied wells. The default `None` draws no labels. |
| `continuous` | bool | False | Use a continuous colormap/colorbar for numeric metrics or gene expression. The default is `False`. |
| `log_values` | bool | True | Apply `log1p` before coloring continuous values. Defaults to `True`; set it to `False` for original or signed scales. |
| `show_legend` | bool | True | Show a categorical legend or continuous colorbar. Categorical entries always represent `plate_value`; fixed-color plots have neither. |
| `nrow` | int | 8 | Number of plate rows. Use 8 for 96-well or 16 for 384-well plates. |
| `ncol` | int | 12 | Number of plate columns. Use 12 for 96-well or 24 for 384-well plates. |
| `nrow_plate` | int \| None | None | Maximum number of plate panels in one figure row. Set `1` to stack plates vertically; `None` arranges panels adaptively. |
| `output_name` | str \| None | None | Filename stem. `None` derives a name from the selected value source. |
| `figsize` | tuple[float, float] \| None | None | Optional whole-figure size in inches; `None` adapts to plate shape and plate count. |
| `well_size` | float | 0.76 | Well diameter relative to one plate-grid cell. |
| `well_border` | bool | True | Whether all wells have visible borders; defaults to `True`. |
| `well_border_color` | str | '#111827' | Matplotlib border color applied when `well_border=True`. |
| `well_border_width` | float | 0.8 | Well-border line width in points. |
| `palette` | PaletteSpec | None | Categorical colors as a Seaborn name, color sequence, category-to-color mapping, or `None`. |
| `cmap` | str | 'RdYlBu\_r' | Continuous Matplotlib colormap. The default `RdYlBu_r` provides a bright blue-to-red scale. |
| `label_color` | str | '#111827' | Color of labels drawn inside occupied wells. |
| `label_size` | float \| None | None | Optional label font size; `None` adapts automatically. |
| `label_style` | str | 'normal' | Label font style, such as `'normal'` or `'italic'`. |
| `label_weight` | str \| int | 'normal' | Label font weight, such as `'normal'` or `'bold'`. |

**Returns / side effects:** Paths for one figure in every requested format. Also writes `<output_name>_plot_table.csv` and `plate_layout_coordinates.csv`. Gene maps additionally write gene-resolution and per-gene-expression tables. Continuous tables retain the original `value` and transformed `plot_value`.

There is no `legend_col`: categorical legend labels and colors both come from `plate_value`. Call this function once per desired figure.

<a id="qc-figures"></a>

## QC report figures

Generate distribution, gene-summary, STARsolo, and complete QC reports.

| Function | Description |
| --- | --- |
| [`drugseqpy.generate_qc_metric_violin_figure()`](#api-generate-qc-metric-violin-figure) | Generate selectable essential QC-metric violin panels, grouped by configurable sample metadata. |
| [`drugseqpy.generate_outlier_score_scatter_figure()`](#api-generate-outlier-score-scatter-figure) | Plot one MAD outlier score per sample, colored by a configurable metadata group. |
| [`drugseqpy.generate_umi_density_figure()`](#api-generate-umi-density-figure) | Generate the per-well total-UMI density figure, grouped by configurable sample metadata. |
| [`drugseqpy.generate_detected_gene_density_figure()`](#api-generate-detected-gene-density-figure) | Generate the per-well detected-gene density figure, grouped by configurable sample metadata. |
| [`drugseqpy.generate_qc_distribution_figures()`](#api-generate-qc-distribution-figures) | Compatibility wrapper that generates the violin panel, outlier scatter, and both sample-level density figures. |
| [`drugseqpy.generate_gene_umi_distribution_figures()`](#api-generate-gene-umi-distribution-figures) | Generate gene-level total-UMI distributions for selected identifier types. |
| [`drugseqpy.generate_marker_gene_figure()`](#api-generate-marker-gene-figure) | Generate an expression panel for an explicit set of marker genes. |
| [`drugseqpy.generate_marker_gene_dotplot_figure()`](#api-generate-marker-gene-dotplot-figure) | Generate a marker dot plot combining expression magnitude and detection prevalence. |
| [`drugseqpy.generate_marker_gene_violin_figure()`](#api-generate-marker-gene-violin-figure) | Generate grouped violin distributions for per-well marker expression. |
| [`drugseqpy.generate_gene_set_violin_figure()`](#api-generate-gene-set-violin-figure) | Generate grouped violin distributions for aggregated per-well gene-set expression. |
| [`drugseqpy.generate_gene_summary_figures()`](#api-generate-gene-summary-figures) | Compatibility wrapper for the focused UMI-distribution and marker-panel functions. |
| [`drugseqpy.generate_star_read_qc_figure()`](#api-generate-star-read-qc-figure) | Generate the optional STARsolo multilayered read-summary figure. |
| [`drugseqpy.generate_star_read_count_violin_figure()`](#api-generate-star-read-count-violin-figure) | Plot matched STAR reads for every analyzed well and annotate group means and medians. |
| [`drugseqpy.generate_star_saturation_figures()`](#api-generate-star-saturation-figures) | Generate optional STARsolo read-count versus detected-gene saturation figures. |
| [`drugseqpy.generate_star_qc_figures()`](#api-generate-star-qc-figures) | Compatibility wrapper for the STARsolo doughnut, per-well read violin, and saturation functions. |
| [`drugseqpy.generate_qc_report()`](#api-generate-qc-report) | Convenience orchestrator that runs all public QC metric and figure groups for batch reporting. |

<a id="qc-figures-details"></a>

### Detailed reference

<a id="api-generate-qc-metric-violin-figure"></a>

### `drugseqpy.generate_qc_metric_violin_figure()`

Generate selectable essential QC-metric violin panels, grouped by configurable sample metadata.

```python
generate_qc_metric_violin_figure(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), report_label: 'str | None' = None, metrics: 'Sequence[str] | None' = None, groupby: 'str' = 'plate_id', split_by: 'str | None' = None, nrow_plate: 'int | None' = None, palette: 'PaletteSpec' = None) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData input. Required sample QC metrics are computed automatically when absent. |
| `output_dir` | str \| Path | required | Directory receiving the violin-panel figure. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `report_label` | str \| None | None | Optional filename prefix. `None` derives a stable prefix from observed plate IDs. |
| `metrics` | Sequence[str] \| None | None | Ordered list of numeric `obs` metrics to plot. `None` uses all four defaults: `mt.pct`, `ribo.pct`, `total_umi`, and `n_genes_det`. |
| `groupby` | str | 'plate\_id' | Column in `dsd.obs` used on the x-axis and for median calculation. Defaults to `plate_id`. |
| `split_by` | str \| None | None | Optional metadata column creating separate figure panels. `None` keeps all samples in the same metric panels. |
| `nrow_plate` | int \| None | None | Maximum number of panels per figure row when faceting; `None` uses the metric-aware layout. |
| `palette` | PaletteSpec | None | Group colors as a Seaborn palette name, color sequence, group-to-color mapping, or `None` for adaptive defaults. |

**Returns / side effects:** Paths for `*_QC_Boxplots` in every requested format. Despite the legacy filename, the panels are violin plots with separate y-axes, individual wells, and median labels placed above the violin bodies.

<a id="api-generate-outlier-score-scatter-figure"></a>

### `drugseqpy.generate_outlier_score_scatter_figure()`

Plot one MAD outlier score per sample, colored by a configurable metadata group.

```python
generate_outlier_score_scatter_figure(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), report_label: 'str | None' = None, groupby: 'str' = 'plate_id', split_by: 'str | None' = None, nrow_plate: 'int | None' = None, palette: 'PaletteSpec' = None, threshold: 'float | None' = None, label_outliers: 'bool' = False, label_col: 'str' = 'well_id', max_labels: 'int' = 10, output_name: 'str | None' = None) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData input. `outlier_score` is computed automatically when absent. |
| `output_dir` | str \| Path | required | Directory receiving the scatter figure and `outlier_score_plot_table.csv`. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `report_label` | str \| None | None | Optional filename prefix. `None` derives a stable prefix from observed plate IDs. |
| `groupby` | str | 'plate\_id' | Column in `dsd.obs` controlling point colors and sample ordering. Defaults to `plate_id`. |
| `split_by` | str \| None | None | Optional metadata column creating separate scatter panels. `None` draws one combined scatter. |
| `nrow_plate` | int \| None | None | Maximum number of split panels per figure row. |
| `palette` | PaletteSpec | None | Point colors as a Seaborn palette name, color sequence, group-to-color mapping, or `None`. |
| `threshold` | float \| None | None | Optional horizontal review threshold. `None` draws no cutoff and does not classify samples. |
| `label_outliers` | bool | False | When enabled, label the highest-scoring samples above `threshold`, or the overall highest scores when no threshold is supplied. |
| `label_col` | str | 'well\_id' | Metadata column used for point labels; falls back to `sample_id` when unavailable. |
| `max_labels` | int | 10 | Maximum number of labels, limiting clutter for large experiments. |
| `output_name` | str \| None | None | Optional custom figure stem. When supplied, it also creates a matching `*_plot_table.csv`, allowing multiple differently grouped outlier plots to coexist. |

**Returns / side effects:** Paths for `*_Outlier_Score_Scatter` in every requested format. Also writes a sample-level plotting table with group, score, plot order, and threshold status.

The MAD score prioritizes unusual multivariate QC profiles. A high value is evidence for review, not an automatic removal decision.

<a id="api-generate-umi-density-figure"></a>

### `drugseqpy.generate_umi_density_figure()`

Generate the per-well total-UMI density figure, grouped by configurable sample metadata.

```python
generate_umi_density_figure(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), report_label: 'str | None' = None, groupby: 'str' = 'plate_id', split_by: 'str | None' = None, nrow_plate: 'int | None' = None, palette: 'PaletteSpec' = None) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData input. `total_umi` is computed automatically when absent. |
| `output_dir` | str \| Path | required | Directory receiving the UMI-density figure. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `report_label` | str \| None | None | Optional filename prefix. `None` derives a stable prefix from observed plate IDs. |
| `groupby` | str | 'plate\_id' | Column in `dsd.obs` defining the overlaid density curves. Defaults to `plate_id`. |
| `split_by` | str \| None | None | Optional metadata column creating separate density panels. `None` overlays all groups in one axes. |
| `nrow_plate` | int \| None | None | Maximum number of split panels per figure row. |
| `palette` | PaletteSpec | None | Curve colors as a Seaborn palette name, color sequence, group-to-color mapping, or `None`. |

**Returns / side effects:** Paths for `*_UMI_Density` in every requested format.

<a id="api-generate-detected-gene-density-figure"></a>

### `drugseqpy.generate_detected_gene_density_figure()`

Generate the per-well detected-gene density figure, grouped by configurable sample metadata.

```python
generate_detected_gene_density_figure(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), report_label: 'str | None' = None, groupby: 'str' = 'plate_id', split_by: 'str | None' = None, nrow_plate: 'int | None' = None, palette: 'PaletteSpec' = None) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData input. `n_genes_det` is computed automatically when absent. |
| `output_dir` | str \| Path | required | Directory receiving the detected-gene-density figure. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `report_label` | str \| None | None | Optional filename prefix. `None` derives a stable prefix from observed plate IDs. |
| `groupby` | str | 'plate\_id' | Column in `dsd.obs` defining the overlaid density curves. Defaults to `plate_id`. |
| `split_by` | str \| None | None | Optional metadata column creating separate density panels. `None` overlays all groups in one axes. |
| `nrow_plate` | int \| None | None | Maximum number of split panels per figure row. |
| `palette` | PaletteSpec | None | Curve colors as a Seaborn palette name, color sequence, group-to-color mapping, or `None`. |

**Returns / side effects:** Paths for `*_Gene_Density` in every requested format.

<a id="api-generate-qc-distribution-figures"></a>

### `drugseqpy.generate_qc_distribution_figures()`

Compatibility wrapper that generates the violin panel, outlier scatter, and both sample-level density figures.

```python
generate_qc_distribution_figures(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), report_label: 'str | None' = None, metrics: 'Sequence[str] | None' = None, groupby: 'str' = 'plate_id', split_by: 'str | None' = None, nrow_plate: 'int | None' = None, palette: 'PaletteSpec' = None, outlier_threshold: 'float | None' = None, label_outliers: 'bool' = False, outlier_label_col: 'str' = 'well_id', max_outlier_labels: 'int' = 10) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData input. Required sample QC metrics are computed automatically when absent. |
| `output_dir` | str \| Path | required | Directory receiving the four figure families and outlier plotting table. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `report_label` | str \| None | None | Optional filename prefix. `None` derives a stable prefix from observed plate IDs. |
| `metrics` | Sequence[str] \| None | None | Optional violin-metric subset forwarded to `generate_qc_metric_violin_figure()`. |
| `groupby` | str | 'plate\_id' | Grouping column forwarded to all four focused visualization functions. |
| `split_by` | str \| None | None | Optional faceting column forwarded to all four focused visualization functions. |
| `nrow_plate` | int \| None | None | Maximum number of split panels per figure row. |
| `palette` | PaletteSpec | None | Categorical colors forwarded to all four focused visualization functions. |
| `outlier_threshold` | float \| None | None | Optional review threshold for the included outlier-score scatter. |
| `label_outliers` | bool | False | Whether to label high-scoring points in the included outlier scatter. |
| `outlier_label_col` | str | 'well\_id' | Metadata column used for outlier point labels. |
| `max_outlier_labels` | int | 10 | Maximum number of outlier labels. |

**Returns / side effects:** Combined paths returned by `generate_qc_metric_violin_figure()`, `generate_outlier_score_scatter_figure()`, `generate_umi_density_figure()`, and `generate_detected_gene_density_figure()`.

Existing scripts can continue to use this wrapper. Tutorials can call the four focused functions so each code block corresponds to one figure.

<a id="api-generate-gene-umi-distribution-figures"></a>

### `drugseqpy.generate_gene_umi_distribution_figures()`

Generate gene-level total-UMI distributions for selected identifier types.

```python
generate_gene_umi_distribution_figures(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), identifier_types: 'str | Sequence[str] | None' = None, report_label: 'str | None' = None, groupby: 'str' = 'plate_id', split_by: 'str | None' = None, nrow_plate: 'int | None' = None, palette: 'PaletteSpec' = None) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData containing raw counts and optional `gene_symbol` annotation. |
| `output_dir` | str \| Path | required | Directory receiving UMI-distribution figures and gene-level source tables. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `identifier_types` | str \| Sequence[str] \| None | None | Identifier distribution(s) to draw: `'ensembl_id'`, `'gene_symbol'`, or a sequence containing both. `None` draws both by default. |
| `report_label` | str \| None | None | Optional filename prefix, otherwise derived from plate IDs. |
| `groupby` | str | 'plate\_id' | Metadata column used to calculate separate per-group gene totals and colored density curves. Defaults to `plate_id`. |
| `split_by` | str \| None | None | Optional metadata column creating separate density panels. `None` overlays the `groupby` curves. |
| `nrow_plate` | int \| None | None | Maximum number of split panels per figure row. |
| `palette` | PaletteSpec | None | Colors for `groupby` density curves as a Seaborn name, color sequence, category mapping, or `None`. |

**Returns / side effects:** Paths for the requested Ensembl-ID and/or gene-symbol distributions. Writes overall gene summaries and `gene_umi_distribution_plot_table.csv` with group-specific totals.

<a id="api-generate-marker-gene-figure"></a>

### `drugseqpy.generate_marker_gene_figure()`

Generate an expression panel for an explicit set of marker genes.

```python
generate_marker_gene_figure(dsd: 'DrugSeqData', output_dir: 'str | Path', *, marker_genes: 'Sequence[str] | None' = ('VIM', 'FN1', 'CDH2', 'CD44', 'EGFR', 'MMP9', 'TGFB1', 'SNAI1', 'TWIST1', 'ITGB1'), formats=('pdf', 'png'), report_label: 'str | None' = None, groupby: 'str' = 'plate_id', split_by: 'str | None' = None, nrow_plate: 'int | None' = None, palette: 'PaletteSpec' = None) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData containing raw counts, plate metadata, compound labels, and gene-symbol annotation. |
| `output_dir` | str \| Path | required | Directory receiving the marker panel and marker source tables. |
| `marker_genes` | Sequence[str] \| None | ('VIM', 'FN1', 'CDH2', 'CD44', 'EGFR', 'MMP9', 'TGFB1', 'SNAI1', 'TWIST1', 'ITGB1') | Ordered marker symbols. The default is the explicit package MDA-MB-231 panel: VIM, FN1, CDH2, CD44, EGFR, MMP9, TGFB1, SNAI1, TWIST1, and ITGB1. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `report_label` | str \| None | None | Optional filename prefix, otherwise derived from plate IDs. |
| `groupby` | str | 'plate\_id' | Metadata column used on the x-axis and for grouped means. Defaults to `plate_id`. |
| `split_by` | str \| None | None | Optional metadata column creating separate gene panels for each split value. |
| `nrow_plate` | int \| None | None | Maximum number of gene/split panels per figure row. |
| `palette` | PaletteSpec | None | Group colors as a Seaborn name, color sequence, category mapping, or `None`. |

**Returns / side effects:** Paths for the marker-expression panel when at least one marker is available. Writes requested-gene availability, sample-level plotting values, and grouped mean-expression tables.

<a id="api-generate-marker-gene-dotplot-figure"></a>

### `drugseqpy.generate_marker_gene_dotplot_figure()`

Generate a marker dot plot combining expression magnitude and detection prevalence.

```python
generate_marker_gene_dotplot_figure(dsd: 'DrugSeqData', output_dir: 'str | Path', *, marker_genes: 'Sequence[str] | None' = ('VIM', 'FN1', 'CDH2', 'CD44', 'EGFR', 'MMP9', 'TGFB1', 'SNAI1', 'TWIST1', 'ITGB1'), formats=('pdf', 'png'), report_label: 'str | None' = None, groupby: 'str' = 'plate_id', split_by: 'str | None' = None, nrow_plate: 'int | None' = None, cmap: 'str' = 'RdYlBu_r', dot_min_size: 'float' = 18, dot_max_size: 'float' = 280) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData containing raw counts, grouping metadata, and gene-symbol annotation. |
| `output_dir` | str \| Path | required | Directory receiving the marker dot plot and its source table. |
| `marker_genes` | Sequence[str] \| None | ('VIM', 'FN1', 'CDH2', 'CD44', 'EGFR', 'MMP9', 'TGFB1', 'SNAI1', 'TWIST1', 'ITGB1') | Ordered marker symbols. The default is the explicit package MDA-MB-231 marker panel. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `report_label` | str \| None | None | Optional filename prefix, otherwise derived from plate IDs. |
| `groupby` | str | 'plate\_id' | Metadata column displayed on the y-axis. Dot statistics are calculated independently for each group. |
| `split_by` | str \| None | None | Optional metadata column creating separate dot-plot panels. |
| `nrow_plate` | int \| None | None | Maximum number of split panels per figure row. |
| `cmap` | str | 'RdYlBu\_r' | Continuous Matplotlib colormap for mean expression. |
| `dot_min_size` | float | 18 | Marker area used for zero-percent detection; must be positive. |
| `dot_max_size` | float | 280 | Marker area used for 100-percent detection; must exceed `dot_min_size`. |

**Returns / side effects:** Paths for the marker dot plot. Writes `marker_dotplot_table.csv` with mean expression, percent-expressing wells, and sample count for every gene/group/facet combination.

<a id="api-generate-marker-gene-violin-figure"></a>

### `drugseqpy.generate_marker_gene_violin_figure()`

Generate grouped violin distributions for per-well marker expression.

```python
generate_marker_gene_violin_figure(dsd: 'DrugSeqData', output_dir: 'str | Path', *, marker_genes: 'Sequence[str] | None' = ('VIM', 'FN1', 'CDH2', 'CD44', 'EGFR', 'MMP9', 'TGFB1', 'SNAI1', 'TWIST1', 'ITGB1'), formats=('pdf', 'png'), report_label: 'str | None' = None, groupby: 'str' = 'plate_id', split_by: 'str | None' = None, nrow_plate: 'int | None' = None, palette: 'PaletteSpec' = None, log_values: 'bool' = False) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData containing raw counts, grouping metadata, and gene-symbol annotation. |
| `output_dir` | str \| Path | required | Directory receiving the marker violin plot and its sample-level source table. |
| `marker_genes` | Sequence[str] \| None | ('VIM', 'FN1', 'CDH2', 'CD44', 'EGFR', 'MMP9', 'TGFB1', 'SNAI1', 'TWIST1', 'ITGB1') | Ordered marker symbols. The default is the explicit package MDA-MB-231 marker panel. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `report_label` | str \| None | None | Optional filename prefix, otherwise derived from plate IDs. |
| `groupby` | str | 'plate\_id' | Metadata column represented by side-by-side colored violins. |
| `split_by` | str \| None | None | Optional metadata column creating a separate gene panel for every split value. |
| `nrow_plate` | int \| None | None | Maximum number of gene/split violin panels per figure row; defaults to five. |
| `palette` | PaletteSpec | None | Group colors as a Seaborn name, color sequence, category mapping, or `None`. |
| `log_values` | bool | False | Apply `log1p` to per-well expression before plotting; defaults to raw expression. |

**Returns / side effects:** Paths for the marker violin plot. Writes `marker_violin_plot_table.csv` with raw and plotted expression for every well and gene.

<a id="api-generate-gene-set-violin-figure"></a>

### `drugseqpy.generate_gene_set_violin_figure()`

Generate grouped violin distributions for aggregated per-well gene-set expression.

```python
generate_gene_set_violin_figure(dsd: 'DrugSeqData', output_dir: 'str | Path', gene_sets: 'Mapping[str, Sequence[str]]', *, gene_aggregation: 'str' = 'mean', formats=('pdf', 'png'), groupby: 'str' = 'plate_id', split_by: 'str | None' = None, nrow_plate: 'int | None' = None, palette: 'PaletteSpec' = None, log_values: 'bool' = True, point_color: 'str' = '#253746', annotation_color: 'str' = '#D32F2F', output_name: 'str' = 'gene_set_expression_violin') -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData containing raw counts, grouping metadata, and gene-symbol annotation. |
| `output_dir` | str \| Path | required | Directory receiving the gene-set violin figure and source tables. |
| `gene_sets` | Mapping[str, Sequence[str]] | required | Mapping from display names to sequences of exact gene symbols or feature IDs. |
| `gene_aggregation` | str | 'mean' | Per-well aggregation across each resolved gene set: `'mean'`, `'sum'`, or `'median'`. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `groupby` | str | 'plate\_id' | Metadata column represented by side-by-side colored violins. |
| `split_by` | str \| None | None | Optional metadata column creating separate gene-set panels. |
| `nrow_plate` | int \| None | None | Maximum number of panels per figure row. |
| `palette` | PaletteSpec | None | Group colors as a Seaborn name, sequence, category mapping, or `None`. |
| `log_values` | bool | True | Apply `log1p` to aggregated expression for plotting while retaining raw values in the source table. |
| `point_color` | str | '#253746' | Color of individual well points overlaid on the violins. |
| `annotation_color` | str | '#D32F2F' | Color of median labels placed above each distribution. |
| `output_name` | str | 'gene\_set\_expression\_violin' | Filename stem shared by the figure and its source tables. |

**Returns / side effects:** Figure paths plus per-well plotting values, gene-resolution records, and per-gene expression tables.

This function uses the same exact identifier resolution and aggregation implementation as `generate_plate_layout_figures()`.

<a id="api-generate-gene-summary-figures"></a>

### `drugseqpy.generate_gene_summary_figures()`

Compatibility wrapper for the focused UMI-distribution and marker-panel functions.

```python
generate_gene_summary_figures(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), identifier_types: 'str | Sequence[str] | None' = None, marker_genes: 'Sequence[str] | None' = ('VIM', 'FN1', 'CDH2', 'CD44', 'EGFR', 'MMP9', 'TGFB1', 'SNAI1', 'TWIST1', 'ITGB1'), report_label: 'str | None' = None, groupby: 'str' = 'plate_id', split_by: 'str | None' = None, nrow_plate: 'int | None' = None, palette: 'PaletteSpec' = None) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData containing raw counts and optional gene-symbol annotation. |
| `output_dir` | str \| Path | required | Directory receiving both requested gene-level figure families. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `identifier_types` | str \| Sequence[str] \| None | None | Identifier distribution(s) forwarded to `generate_gene_umi_distribution_figures()`; `None` selects both. |
| `marker_genes` | Sequence[str] \| None | ('VIM', 'FN1', 'CDH2', 'CD44', 'EGFR', 'MMP9', 'TGFB1', 'SNAI1', 'TWIST1', 'ITGB1') | Markers forwarded to `generate_marker_gene_figure()`. The explicit MDA-MB-231 panel is the default; an empty sequence disables this figure. |
| `report_label` | str \| None | None | Optional filename prefix, otherwise derived from plate IDs. |
| `groupby` | str | 'plate\_id' | Metadata column forwarded to both focused functions; defaults to `plate_id`. |
| `split_by` | str \| None | None | Optional faceting column forwarded to both focused functions. |
| `nrow_plate` | int \| None | None | Maximum number of split panels per figure row. |
| `palette` | PaletteSpec | None | Palette forwarded to both focused functions. |

**Returns / side effects:** Combined figure paths returned by the two focused functions.

New tutorials should call the focused functions separately so each call maps to one figure family.

<a id="api-generate-star-read-qc-figure"></a>

### `drugseqpy.generate_star_read_qc_figure()`

Generate the optional STARsolo multilayered read-summary figure.

```python
generate_star_read_qc_figure(dsd: 'DrugSeqData', star_root: 'str | Path', output_dir: 'str | Path', *, formats=('pdf', 'png'), report_label: 'str | None' = None, groupby: 'str' = 'plate_id', split_by: 'str | None' = None, nrow_plate: 'int | None' = None, palette: 'PaletteSpec' = None) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData used to identify plates and their grouping values. |
| `star_root` | str \| Path | required | Generic STAR root containing `<plate>/library_seq_Solo.out/Gene/Summary.csv`. |
| `output_dir` | str \| Path | required | Directory receiving the doughnut figure and `star_group_summary.csv`. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `report_label` | str \| None | None | Optional filename prefix, otherwise derived from plate IDs. |
| `groupby` | str | 'plate\_id' | Metadata column used to label/group panels. It must resolve to one value per plate and defaults to `plate_id`. |
| `split_by` | str \| None | None | Optional plate-level metadata column used instead of `groupby` to create separate doughnut panels. |
| `nrow_plate` | int \| None | None | Maximum number of grouped plate panels per figure row. `1` stacks panels vertically. |
| `palette` | PaletteSpec | None | Colors for doughnut rings as a Seaborn name, color sequence, mapping, or `None`. |

**Returns / side effects:** Paths for the multilayered doughnut figure. Missing `Summary.csv` inputs are skipped rather than fabricated.

<a id="api-generate-star-read-count-violin-figure"></a>

### `drugseqpy.generate_star_read_count_violin_figure()`

Plot matched STAR reads for every analyzed well and annotate group means and medians.

```python
generate_star_read_count_violin_figure(dsd: 'DrugSeqData', star_root: 'str | Path', output_dir: 'str | Path', *, formats=('pdf', 'png'), report_label: 'str | None' = None, groupby: 'str' = 'plate_id', split_by: 'str | None' = None, nrow_plate: 'int | None' = None, palette: 'PaletteSpec' = None, point_color: 'str' = '#253746', annotation_color: 'str' = '#D32F2F') -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData defining the full analyzed-well set and plate/group metadata. |
| `star_root` | str \| Path | required | STAR root containing per-plate `CellReads.stats` and filtered barcodes. |
| `output_dir` | str \| Path | required | Directory receiving the violin figure, per-well read table, and group summary table. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `report_label` | str \| None | None | Optional filename prefix, otherwise derived from observed plate IDs. |
| `groupby` | str | 'plate\_id' | Metadata column used for violin fills, x-axis groups, means, and medians; defaults to `plate_id`. |
| `split_by` | str \| None | None | Optional metadata column creating separate violin panels. |
| `nrow_plate` | int \| None | None | Maximum number of split panels per figure row. |
| `palette` | PaletteSpec | None | Violin fill colors as a Seaborn name, color sequence, category mapping, or `None`. |
| `point_color` | str | '#253746' | Color of individual well points overlaid on the violins. |
| `annotation_color` | str | '#D32F2F' | Color of mean and median labels above each violin. |

**Returns / side effects:** Paths for `*_STAR_Read_Count_Violin`. Writes `star_read_count_per_well.csv` and `star_read_count_summary.csv`; unmatched wells are retained with zero reads and an explicit match flag.

<a id="api-generate-star-saturation-figures"></a>

### `drugseqpy.generate_star_saturation_figures()`

Generate optional STARsolo read-count versus detected-gene saturation figures.

```python
generate_star_saturation_figures(dsd: 'DrugSeqData', star_root: 'str | Path', output_dir: 'str | Path', *, formats=('pdf', 'png'), report_label: 'str | None' = None, groupby: 'str' = 'plate_id', split_by: 'str | None' = None, nrow_plate: 'int | None' = None, palette: 'PaletteSpec' = None) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData used to connect STAR barcodes to matrix-level detected-gene counts. |
| `star_root` | str \| Path | required | Generic STAR root containing per-plate `CellReads.stats` and, preferably, filtered barcodes. |
| `output_dir` | str \| Path | required | Directory receiving saturation figures and `star_saturation_metrics.csv`. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `report_label` | str \| None | None | Optional filename prefix, otherwise derived from plate IDs. |
| `groupby` | str | 'plate\_id' | Metadata column defining point colors in each saturation axes; defaults to `plate_id`. |
| `split_by` | str \| None | None | Optional metadata column creating separate saturation panels. `None` combines all samples in one axes. |
| `nrow_plate` | int \| None | None | Maximum number of split saturation panels per figure row. `1` stacks panels vertically. |
| `palette` | PaletteSpec | None | Panel colors as a Seaborn name, color sequence, mapping, or `None`. |

**Returns / side effects:** Paths for available Ensembl-ID and gene-symbol saturation plots. Axes use original, non-log read and detected-gene values.

<a id="api-generate-star-qc-figures"></a>

### `drugseqpy.generate_star_qc_figures()`

Compatibility wrapper for the STARsolo doughnut, per-well read violin, and saturation functions.

```python
generate_star_qc_figures(dsd: 'DrugSeqData', star_root: 'str | Path', output_dir: 'str | Path', *, formats=('pdf', 'png'), report_label: 'str | None' = None, groupby: 'str' = 'plate_id', split_by: 'str | None' = None, nrow_plate: 'int | None' = None, palette: 'PaletteSpec' = None) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData used to identify plate order and connect barcodes to matrix-level detected-gene counts. |
| `star_root` | str \| Path | required | STAR root containing `<plate>/library_seq_Solo.out/Gene` directories. |
| `output_dir` | str \| Path | required | Directory receiving STAR figures and source tables. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. |
| `report_label` | str \| None | None | Optional filename prefix, otherwise derived from plate IDs. |
| `groupby` | str | 'plate\_id' | Metadata column defining figure groups; defaults to `plate_id`. |
| `split_by` | str \| None | None | Optional metadata column defining separate panels. |
| `nrow_plate` | int \| None | None | Maximum number of grouped panels per figure row. |
| `palette` | PaletteSpec | None | Colors for doughnut rings and plate groups as a Seaborn name, color sequence, category mapping, or `None`. |

**Returns / side effects:** Paths for the multilayered doughnut, per-well read-count violin, and available Ensembl-ID/gene-symbol saturation plots. Missing STAR inputs are skipped rather than fabricated.

<a id="api-generate-qc-report"></a>

### `drugseqpy.generate_qc_report()`

Convenience orchestrator that runs all public QC metric and figure groups for batch reporting.

```python
generate_qc_report(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), star_root: 'str | Path | None' = None, marker_genes: 'Sequence[str] | None' = ('VIM', 'FN1', 'CDH2', 'CD44', 'EGFR', 'MMP9', 'TGFB1', 'SNAI1', 'TWIST1', 'ITGB1'), gene_umi_identifier_types: 'str | Sequence[str] | None' = None, plate_gene: 'str | Sequence[str] | None' = 'FN1', plate_gene_set: 'Sequence[str] | None' = None, plate_gene_aggregation: 'str' = 'mean', plate_nrow: 'int' = 8, plate_ncol: 'int' = 12, nrow_plate: 'int | None' = None, report_label: 'str | None' = None, violin_metrics: 'Sequence[str] | None' = None, groupby: 'str' = 'plate_id', split_by: 'str | None' = None, palette: 'PaletteSpec' = None, cmap: 'str' = 'RdYlBu_r', outlier_threshold: 'float | None' = None, label_outliers: 'bool' = False, outlier_label_col: 'str' = 'well_id', max_outlier_labels: 'int' = 10) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | Merged DrugSeqData input. |
| `output_dir` | str \| Path | required | Directory receiving figures, CSV source tables, and the figure manifest. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib-supported output extensions. This HTML tutorial requests PNG only. |
| `star_root` | str \| Path \| None | None | Optional STAR result root containing `<plate>/library_seq_Solo.out/Gene`. If absent, STAR-only plots are skipped. |
| `marker_genes` | Sequence[str] \| None | ('VIM', 'FN1', 'CDH2', 'CD44', 'EGFR', 'MMP9', 'TGFB1', 'SNAI1', 'TWIST1', 'ITGB1') | Marker symbols for the focused marker panel. The default is the explicit MDA-MB-231 panel: VIM, FN1, CDH2, CD44, EGFR, MMP9, TGFB1, SNAI1, TWIST1, and ITGB1. |
| `gene_umi_identifier_types` | str \| Sequence[str] \| None | None | UMI-distribution identifier types. `None` draws both Ensembl-ID and gene-symbol distributions; pass one type to select only that figure. |
| `plate_gene` | str \| Sequence[str] \| None | 'FN1' | Optional string gene or gene sequence used for the primary spatial expression map; defaults to `FN1` and is skipped if unavailable. |
| `plate_gene_set` | Sequence[str] \| None | None | Optional additional gene sequence used to generate a separate gene-set expression plate map. |
| `plate_gene_aggregation` | str | 'mean' | Aggregation used when either expression specification is a sequence: `'mean'`, `'sum'`, or `'median'`. |
| `plate_nrow` | int | 8 | Number of rows used for all report plate layouts; defaults to 8. |
| `plate_ncol` | int | 12 | Number of columns used for all report plate layouts; defaults to 12. |
| `nrow_plate` | int \| None | None | Maximum plate/group panels per row for plate layouts and STARsolo figures. |
| `report_label` | str \| None | None | Optional filename prefix. `None` derives a combined label from the observed plate IDs. |
| `violin_metrics` | Sequence[str] \| None | None | Optional ordered subset of QC metrics shown in the violin panel; `None` uses all four defaults. |
| `groupby` | str | 'plate\_id' | Metadata column used by the QC violin, outlier-score scatter, and density figures. Defaults to `plate_id`. |
| `split_by` | str \| None | None | Optional metadata column creating separate panels across group-aware QC visualizations. `None` keeps combined plots. |
| `palette` | PaletteSpec | None | Categorical palette forwarded to all QC visualization functions. Accepts a Seaborn name, color sequence, or category-to-color mapping. |
| `cmap` | str | 'RdYlBu\_r' | Continuous Matplotlib colormap forwarded to numeric plate layouts. |
| `outlier_threshold` | float \| None | None | Optional review threshold forwarded to the outlier-score scatter. `None` draws no threshold line. |
| `label_outliers` | bool | False | Label high-scoring points in the outlier-score scatter when enabled. |
| `outlier_label_col` | str | 'well\_id' | Metadata column used for labels in the outlier-score scatter. |
| `max_outlier_labels` | int | 10 | Maximum number of outlier labels drawn in the complete report. |

**Returns / side effects:** A list of all figure paths. The object gains QC columns; `qc_figure_manifest.csv` records the public producer function for every figure.

Use this wrapper for a complete automated report. In a tutorial or interactive analysis, call the focused figure functions separately so each code block has a visible, inspectable result.

<a id="mix-species"></a>

## Mix-species analysis

Configure any two species, assign wells, and quantify reciprocal cross-contamination.

| Function | Description |
| --- | --- |
| [`drugseqpy.assign_well_species()`](#api-assign-well-species) | Assign the designed base species of each well from an existing metadata column. |
| [`drugseqpy.compute_mix_species_qc()`](#api-compute-mix-species-qc) | Count any two configured species, classify wells, and quantify cross-well contamination relative to the assigned base species. |
| [`drugseqpy.generate_mix_species_barnyard_figure()`](#api-generate-mix-species-barnyard-figure) | Generate only the plate-faceted species-A-versus-species-B UMI barnyard figure with square 1:1 coordinate panels. |
| [`drugseqpy.generate_mix_species_purity_figure()`](#api-generate-mix-species-purity-figure) | Generate only the plate-faceted two-species purity figure. |
| [`drugseqpy.generate_mix_species_concordance_figure()`](#api-generate-mix-species-concordance-figure) | Generate only the expected-label versus observed-species stacked-bar figure. |
| [`drugseqpy.generate_mix_species_qc_metrics_figure()`](#api-generate-mix-species-qc-metrics-figure) | Generate only the plate-level mix-species technical-QC violin panel. |
| [`drugseqpy.generate_mix_species_fraction_figure()`](#api-generate-mix-species-fraction-figure) | Generate only the two-species fraction-by-expected-label figure. |
| [`drugseqpy.generate_mix_species_contamination_figure()`](#api-generate-mix-species-contamination-figure) | Plot off-target fractions separately with species A and species B as the assigned base. |
| [`drugseqpy.generate_mix_species_report()`](#api-generate-mix-species-report) | Generate the independent mix-species figure set and all plot source/audit tables. |

<a id="mix-species-details"></a>

### Detailed reference

<a id="api-assign-well-species"></a>

### `drugseqpy.assign_well_species()`

Assign the designed base species of each well from an existing metadata column.

```python
assign_well_species(dsd: 'DrugSeqData', source_col: 'str', mapping: 'Mapping[str, str]', *, output_col: 'str' = 'expected_species', default: 'str' = 'unknown', inplace: 'bool' = True) -> 'DrugSeqData | None'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData whose `obs` contains the source labels. |
| `source_col` | str | required | Metadata column containing cell-line, sample-class, or other design labels. |
| `mapping` | Mapping[str, str] | required | Mapping from each design label to one configured species name, `control`, or `unknown`. |
| `output_col` | str | 'expected\_species' | Metadata column that receives the assigned base species. |
| `default` | str | 'unknown' | Value assigned to source labels absent from `mapping`. |
| `inplace` | bool | True | Modify the input object when `True`; otherwise return a modified copy. |

**Returns / side effects:** Stores the well assignment in `obs[output_col]` and the mapping provenance in `uns['well_species_assignment']`; returns a copy only when `inplace=False`.

Assignment is case-insensitive. The two species need not be human and mouse.

<a id="api-compute-mix-species-qc"></a>

### `drugseqpy.compute_mix_species_qc()`

Count any two configured species, classify wells, and quantify cross-well contamination relative to the assigned base species.

```python
compute_mix_species_qc(dsd: 'DrugSeqData', *, species_a_prefix: 'str' = 'ENSG', species_b_prefix: 'str' = 'ENSMUSG', species_a_name: 'str' = 'human', species_b_name: 'str' = 'mouse', expected_col: 'str | None' = None, expected_species_map: 'Mapping[str, str] | None' = None, assignment_threshold: 'float' = 0.5, min_species_umi: 'float' = 0, mitochondrial_pattern: 'str' = '^(?:MT-|mt-)', human_prefix: 'str | None' = None, mouse_prefix: 'str | None' = None, human_threshold: 'float | None' = None, min_total_umi: 'float | None' = None) -> 'pd.DataFrame'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData with raw two-species counts. |
| `species_a_prefix` | str | 'ENSG' | Gene-ID prefix selecting species A features. |
| `species_b_prefix` | str | 'ENSMUSG' | Gene-ID prefix selecting species B features. |
| `species_a_name` | str | 'human' | Display and assignment name for species A. |
| `species_b_name` | str | 'mouse' | Display and assignment name for species B. |
| `expected_col` | str \| None | None | Optional metadata column containing the assigned base species, or labels that can be resolved with `expected_species_map`. |
| `expected_species_map` | Mapping[str, str] \| None | None | Optional mapping from values in `expected_col` to species A, species B, `control`, or `unknown`. |
| `assignment_threshold` | float | 0.5 | Minimum within-pair fraction required for a species call. Use a value above 0.5 for mutually exclusive majority criteria; at a 0.5 tie, species A is selected first. This is a classification threshold, not a validated purity cutoff. |
| `min_species_umi` | float | 0 | Minimum species-A-plus-species-B UMI count; lower-count wells are called `low_umi`. |
| `mitochondrial_pattern` | str | '^(?:MT-\|mt-)' | Regular expression matched against gene symbols to calculate mitochondrial UMI percentage. |
| `human_prefix` | str \| None | None | Backward-compatible alias for `species_a_prefix`. |
| `mouse_prefix` | str \| None | None | Backward-compatible alias for `species_b_prefix`. |
| `human_threshold` | float \| None | None | Backward-compatible alias for `assignment_threshold`. |
| `min_total_umi` | float \| None | None | Backward-compatible alias for `min_species_umi`. |

**Returns / side effects:** A per-well DataFrame with generic species-A/B counts and fractions, observed calls, and—when assignments are supplied—base species, on-target signal, off-target signal, cross-contamination fraction, and concordance.

Human–mouse defaults and historical column aliases are retained for compatibility, but every new analysis can provide arbitrary species names and feature prefixes.

<a id="api-generate-mix-species-barnyard-figure"></a>

### `drugseqpy.generate_mix_species_barnyard_figure()`

Generate only the plate-faceted species-A-versus-species-B UMI barnyard figure with square 1:1 coordinate panels.

```python
generate_mix_species_barnyard_figure(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), expected_col: 'str | None' = 'compound', species_a_prefix: 'str' = 'ENSG', species_b_prefix: 'str' = 'ENSMUSG', species_a_name: 'str' = 'human', species_b_name: 'str' = 'mouse', expected_species_map: 'Mapping[str, str] | None' = None, assignment_threshold: 'float' = 0.5, min_species_umi: 'float' = 0, human_prefix: 'str | None' = None, mouse_prefix: 'str | None' = None, nrow_plate: 'int | None' = None, metric_species: 'str | Sequence[str] | None' = 'base', expected_palette: 'PaletteSpec' = None, reference_line_color: 'str' = '#7F7F7F', annotation_color: 'str' = '#D32F2F', edge_color: 'str' = 'white', grid_color: 'str' = '#DCE6F2', panel_border_color: 'str' = '#334155', point_alpha: 'float' = 0.72, show_observed_style: 'bool' = True) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData with raw two-species counts. |
| `output_dir` | str \| Path | required | Directory receiving the figure and its plot-source table. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib output formats. |
| `expected_col` | str \| None | 'compound' | Optional metadata column containing the assigned base species, or labels that can be resolved with `expected_species_map`. |
| `species_a_prefix` | str | 'ENSG' | Gene-ID prefix selecting species A features. |
| `species_b_prefix` | str | 'ENSMUSG' | Gene-ID prefix selecting species B features. |
| `species_a_name` | str | 'human' | Display and assignment name for species A. |
| `species_b_name` | str | 'mouse' | Display and assignment name for species B. |
| `expected_species_map` | Mapping[str, str] \| None | None | Optional mapping from values in `expected_col` to species A, species B, `control`, or `unknown`. |
| `assignment_threshold` | float | 0.5 | Minimum within-pair fraction required for a species call. Use a value above 0.5 for mutually exclusive majority criteria; at a 0.5 tie, species A is selected first. This is a classification threshold, not a validated purity cutoff. |
| `min_species_umi` | float | 0 | Minimum species-A-plus-species-B UMI count; lower-count wells are called `low_umi`. |
| `human_prefix` | str \| None | None | Backward-compatible alias for `species_a_prefix`. |
| `mouse_prefix` | str \| None | None | Backward-compatible alias for `species_b_prefix`. |
| `nrow_plate` | int \| None | None | Maximum number of plate facets per figure row. |
| `metric_species` | str \| Sequence[str] \| None | 'base' | Base-species purity summaries to annotate: `'base'` (default, all assigned bases present), `'both'`, `'species_a'`, `'species_b'`, one configured species name, a sequence of names, or `None` to hide metrics. |
| `expected_palette` | PaletteSpec | None | Colors for expected labels or cell types. |
| `reference_line_color` | str | '#7F7F7F' | Color of the equal-species-count diagonal. |
| `annotation_color` | str | '#D32F2F' | Color of per-panel summary annotations. |
| `edge_color` | str | 'white' | Point-edge color. |
| `grid_color` | str | '#DCE6F2' | Dashed panel-grid color. |
| `panel_border_color` | str | '#334155' | Complete panel-border color. |
| `point_alpha` | float | 0.72 | Point opacity between zero and one. |
| `show_observed_style` | bool | True | Whether point marker shape also represents observed species. |

**Returns / side effects:** Paths for `species_barnyard` in every requested format; also writes the per-well plot table and `species_barnyard_metric_summary.csv`.

Each facet uses identical x/y limits and an enforced square plotting box. Selected base-species annotations report minimum, mean, median, and maximum purity.

<a id="api-generate-mix-species-purity-figure"></a>

### `drugseqpy.generate_mix_species_purity_figure()`

Generate only the plate-faceted two-species purity figure.

```python
generate_mix_species_purity_figure(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), expected_col: 'str | None' = 'compound', species_a_prefix: 'str' = 'ENSG', species_b_prefix: 'str' = 'ENSMUSG', species_a_name: 'str' = 'human', species_b_name: 'str' = 'mouse', expected_species_map: 'Mapping[str, str] | None' = None, assignment_threshold: 'float' = 0.5, min_species_umi: 'float' = 0, human_prefix: 'str | None' = None, mouse_prefix: 'str | None' = None, nrow_plate: 'int | None' = None, expected_palette: 'PaletteSpec' = None, reference_line_color: 'str' = '#7F7F7F', guide_line_color: 'str' = '#CBD5E1', edge_color: 'str' = 'white', grid_color: 'str' = '#DCE6F2', panel_border_color: 'str' = '#334155', point_alpha: 'float' = 0.72, show_observed_style: 'bool' = True) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData with raw two-species counts. |
| `output_dir` | str \| Path | required | Directory receiving the figure and its plot-source table. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib output formats. |
| `expected_col` | str \| None | 'compound' | Optional metadata column containing the assigned base species, or labels that can be resolved with `expected_species_map`. |
| `species_a_prefix` | str | 'ENSG' | Gene-ID prefix selecting species A features. |
| `species_b_prefix` | str | 'ENSMUSG' | Gene-ID prefix selecting species B features. |
| `species_a_name` | str | 'human' | Display and assignment name for species A. |
| `species_b_name` | str | 'mouse' | Display and assignment name for species B. |
| `expected_species_map` | Mapping[str, str] \| None | None | Optional mapping from values in `expected_col` to species A, species B, `control`, or `unknown`. |
| `assignment_threshold` | float | 0.5 | Minimum within-pair fraction required for a species call. Use a value above 0.5 for mutually exclusive majority criteria; at a 0.5 tie, species A is selected first. This is a classification threshold, not a validated purity cutoff. |
| `min_species_umi` | float | 0 | Minimum species-A-plus-species-B UMI count; lower-count wells are called `low_umi`. |
| `human_prefix` | str \| None | None | Backward-compatible alias for `species_a_prefix`. |
| `mouse_prefix` | str \| None | None | Backward-compatible alias for `species_b_prefix`. |
| `nrow_plate` | int \| None | None | Maximum number of plate facets per figure row. |
| `expected_palette` | PaletteSpec | None | Colors for expected labels or cell types. |
| `reference_line_color` | str | '#7F7F7F' | Color of the complementary-fraction diagonal. |
| `guide_line_color` | str | '#CBD5E1' | Color of the 50% guides. |
| `edge_color` | str | 'white' | Point-edge color. |
| `grid_color` | str | '#DCE6F2' | Dashed panel-grid color. |
| `panel_border_color` | str | '#334155' | Complete panel-border color. |
| `point_alpha` | float | 0.72 | Point opacity between zero and one. |
| `show_observed_style` | bool | True | Whether point marker shape also represents observed species. |

**Returns / side effects:** Paths for `species_purity` in every requested format; also writes `species_purity_plot_table.csv`.

Every facet has identical 0–100% x/y limits and an enforced square 1:1 plotting box.

<a id="api-generate-mix-species-concordance-figure"></a>

### `drugseqpy.generate_mix_species_concordance_figure()`

Generate only the expected-label versus observed-species stacked-bar figure.

```python
generate_mix_species_concordance_figure(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), expected_col: 'str | None' = 'compound', species_a_prefix: 'str' = 'ENSG', species_b_prefix: 'str' = 'ENSMUSG', species_a_name: 'str' = 'human', species_b_name: 'str' = 'mouse', expected_species_map: 'Mapping[str, str] | None' = None, assignment_threshold: 'float' = 0.5, min_species_umi: 'float' = 0, human_prefix: 'str | None' = None, mouse_prefix: 'str | None' = None, nrow_plate: 'int | None' = None, observed_palette: 'PaletteSpec' = None, edge_color: 'str' = 'white', label_color: 'str' = 'white', grid_color: 'str' = '#DCE6F2', panel_border_color: 'str' = '#334155', fill_alpha: 'float' = 0.95) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData with raw two-species counts. |
| `output_dir` | str \| Path | required | Directory receiving the figure and its plot-source table. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib output formats. |
| `expected_col` | str \| None | 'compound' | Optional metadata column containing the assigned base species, or labels that can be resolved with `expected_species_map`. |
| `species_a_prefix` | str | 'ENSG' | Gene-ID prefix selecting species A features. |
| `species_b_prefix` | str | 'ENSMUSG' | Gene-ID prefix selecting species B features. |
| `species_a_name` | str | 'human' | Display and assignment name for species A. |
| `species_b_name` | str | 'mouse' | Display and assignment name for species B. |
| `expected_species_map` | Mapping[str, str] \| None | None | Optional mapping from values in `expected_col` to species A, species B, `control`, or `unknown`. |
| `assignment_threshold` | float | 0.5 | Minimum within-pair fraction required for a species call. Use a value above 0.5 for mutually exclusive majority criteria; at a 0.5 tie, species A is selected first. This is a classification threshold, not a validated purity cutoff. |
| `min_species_umi` | float | 0 | Minimum species-A-plus-species-B UMI count; lower-count wells are called `low_umi`. |
| `human_prefix` | str \| None | None | Backward-compatible alias for `species_a_prefix`. |
| `mouse_prefix` | str \| None | None | Backward-compatible alias for `species_b_prefix`. |
| `nrow_plate` | int \| None | None | Maximum number of plate facets per figure row. |
| `observed_palette` | PaletteSpec | None | Fill colors for observed-species categories. |
| `edge_color` | str | 'white' | Stacked-bar edge color. |
| `label_color` | str | 'white' | Color of counts printed inside bars. |
| `grid_color` | str | '#DCE6F2' | Dashed panel-grid color. |
| `panel_border_color` | str | '#334155' | Complete panel-border color. |
| `fill_alpha` | float | 0.95 | Bar opacity between zero and one. |

**Returns / side effects:** Paths for `species_concordance_counts`; also writes `species_concordance_plot_table.csv`.

<a id="api-generate-mix-species-qc-metrics-figure"></a>

### `drugseqpy.generate_mix_species_qc_metrics_figure()`

Generate only the plate-level mix-species technical-QC violin panel.

```python
generate_mix_species_qc_metrics_figure(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), expected_col: 'str | None' = 'compound', species_a_prefix: 'str' = 'ENSG', species_b_prefix: 'str' = 'ENSMUSG', species_a_name: 'str' = 'human', species_b_name: 'str' = 'mouse', expected_species_map: 'Mapping[str, str] | None' = None, assignment_threshold: 'float' = 0.5, min_species_umi: 'float' = 0, human_prefix: 'str | None' = None, mouse_prefix: 'str | None' = None, plate_palette: 'PaletteSpec' = None, point_color: 'str' = '#253746', annotation_color: 'str' = '#D32F2F', grid_color: 'str' = '#DCE6F2', panel_border_color: 'str' = '#334155') -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData with raw two-species counts. |
| `output_dir` | str \| Path | required | Directory receiving the figure and its plot-source table. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib output formats. |
| `expected_col` | str \| None | 'compound' | Optional metadata column containing the assigned base species, or labels that can be resolved with `expected_species_map`. |
| `species_a_prefix` | str | 'ENSG' | Gene-ID prefix selecting species A features. |
| `species_b_prefix` | str | 'ENSMUSG' | Gene-ID prefix selecting species B features. |
| `species_a_name` | str | 'human' | Display and assignment name for species A. |
| `species_b_name` | str | 'mouse' | Display and assignment name for species B. |
| `expected_species_map` | Mapping[str, str] \| None | None | Optional mapping from values in `expected_col` to species A, species B, `control`, or `unknown`. |
| `assignment_threshold` | float | 0.5 | Minimum within-pair fraction required for a species call. Use a value above 0.5 for mutually exclusive majority criteria; at a 0.5 tie, species A is selected first. This is a classification threshold, not a validated purity cutoff. |
| `min_species_umi` | float | 0 | Minimum species-A-plus-species-B UMI count; lower-count wells are called `low_umi`. |
| `human_prefix` | str \| None | None | Backward-compatible alias for `species_a_prefix`. |
| `mouse_prefix` | str \| None | None | Backward-compatible alias for `species_b_prefix`. |
| `plate_palette` | PaletteSpec | None | Colors for plate groups. |
| `point_color` | str | '#253746' | Color of individual observations. |
| `annotation_color` | str | '#D32F2F' | Color of median labels. |
| `grid_color` | str | '#DCE6F2' | Dashed panel-grid color. |
| `panel_border_color` | str | '#334155' | Complete panel-border color. |

**Returns / side effects:** Paths for `mix_species_qc_metrics`; also writes `mix_species_qc_metrics_plot_table.csv`.

<a id="api-generate-mix-species-fraction-figure"></a>

### `drugseqpy.generate_mix_species_fraction_figure()`

Generate only the two-species fraction-by-expected-label figure.

```python
generate_mix_species_fraction_figure(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), expected_col: 'str | None' = 'compound', species_a_prefix: 'str' = 'ENSG', species_b_prefix: 'str' = 'ENSMUSG', species_a_name: 'str' = 'human', species_b_name: 'str' = 'mouse', expected_species_map: 'Mapping[str, str] | None' = None, assignment_threshold: 'float' = 0.5, min_species_umi: 'float' = 0, human_prefix: 'str | None' = None, mouse_prefix: 'str | None' = None, plate_palette: 'PaletteSpec' = None, edge_color: 'str' = 'white', grid_color: 'str' = '#DCE6F2', panel_border_color: 'str' = '#334155', fill_alpha: 'float' = 0.95, point_alpha: 'float' = 0.35) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData with raw two-species counts. |
| `output_dir` | str \| Path | required | Directory receiving the figure and its plot-source table. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib output formats. |
| `expected_col` | str \| None | 'compound' | Optional metadata column containing the assigned base species, or labels that can be resolved with `expected_species_map`. |
| `species_a_prefix` | str | 'ENSG' | Gene-ID prefix selecting species A features. |
| `species_b_prefix` | str | 'ENSMUSG' | Gene-ID prefix selecting species B features. |
| `species_a_name` | str | 'human' | Display and assignment name for species A. |
| `species_b_name` | str | 'mouse' | Display and assignment name for species B. |
| `expected_species_map` | Mapping[str, str] \| None | None | Optional mapping from values in `expected_col` to species A, species B, `control`, or `unknown`. |
| `assignment_threshold` | float | 0.5 | Minimum within-pair fraction required for a species call. Use a value above 0.5 for mutually exclusive majority criteria; at a 0.5 tie, species A is selected first. This is a classification threshold, not a validated purity cutoff. |
| `min_species_umi` | float | 0 | Minimum species-A-plus-species-B UMI count; lower-count wells are called `low_umi`. |
| `human_prefix` | str \| None | None | Backward-compatible alias for `species_a_prefix`. |
| `mouse_prefix` | str \| None | None | Backward-compatible alias for `species_b_prefix`. |
| `plate_palette` | PaletteSpec | None | Colors for plate groups. |
| `edge_color` | str | 'white' | Point-edge color. |
| `grid_color` | str | '#DCE6F2' | Dashed panel-grid color. |
| `panel_border_color` | str | '#334155' | Complete panel-border color. |
| `fill_alpha` | float | 0.95 | Boxplot opacity between zero and one. |
| `point_alpha` | float | 0.35 | Individual-point opacity between zero and one. |

**Returns / side effects:** Paths for `species_fraction_by_expected`; also writes `species_fraction_plot_table.csv`.

<a id="api-generate-mix-species-contamination-figure"></a>

### `drugseqpy.generate_mix_species_contamination_figure()`

Plot off-target fractions separately with species A and species B as the assigned base.

```python
generate_mix_species_contamination_figure(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), expected_col: 'str | None' = 'expected_species', species_a_prefix: 'str' = 'ENSG', species_b_prefix: 'str' = 'ENSMUSG', species_a_name: 'str' = 'human', species_b_name: 'str' = 'mouse', expected_species_map: 'Mapping[str, str] | None' = None, assignment_threshold: 'float' = 0.5, min_species_umi: 'float' = 0, human_prefix: 'str | None' = None, mouse_prefix: 'str | None' = None, plate_palette: 'PaletteSpec' = None, point_color: 'str' = '#253746', annotation_color: 'str' = '#D32F2F', grid_color: 'str' = '#DCE6F2', panel_border_color: 'str' = '#334155') -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | DrugSeqData with raw two-species counts and well assignments. |
| `output_dir` | str \| Path | required | Directory receiving the figure and its plot-source table. |
| `formats` | — | ('pdf', 'png') | Iterable of Matplotlib output formats. |
| `expected_col` | str \| None | 'expected\_species' | Optional metadata column containing the assigned base species, or labels that can be resolved with `expected_species_map`. |
| `species_a_prefix` | str | 'ENSG' | Gene-ID prefix selecting species A features. |
| `species_b_prefix` | str | 'ENSMUSG' | Gene-ID prefix selecting species B features. |
| `species_a_name` | str | 'human' | Display and assignment name for species A. |
| `species_b_name` | str | 'mouse' | Display and assignment name for species B. |
| `expected_species_map` | Mapping[str, str] \| None | None | Optional mapping from values in `expected_col` to species A, species B, `control`, or `unknown`. |
| `assignment_threshold` | float | 0.5 | Minimum within-pair fraction required for a species call. Use a value above 0.5 for mutually exclusive majority criteria; at a 0.5 tie, species A is selected first. This is a classification threshold, not a validated purity cutoff. |
| `min_species_umi` | float | 0 | Minimum species-A-plus-species-B UMI count; lower-count wells are called `low_umi`. |
| `human_prefix` | str \| None | None | Backward-compatible alias for `species_a_prefix`. |
| `mouse_prefix` | str \| None | None | Backward-compatible alias for `species_b_prefix`. |
| `plate_palette` | PaletteSpec | None | Colors for plate groups. |
| `point_color` | str | '#253746' | Color of individual wells. |
| `annotation_color` | str | '#D32F2F' | Color of median annotations. |
| `grid_color` | str | '#DCE6F2' | Dashed panel-grid color. |
| `panel_border_color` | str | '#334155' | Complete panel-border color. |

**Returns / side effects:** Paths for `species_cross_contamination_by_base`; also writes its per-well plot table.

The left panel treats species A as the intended species, so species B is off target; the right panel reverses that direction. Both fractions divide by species-A-plus-species-B UMIs.

<a id="api-generate-mix-species-report"></a>

### `drugseqpy.generate_mix_species_report()`

Generate the independent mix-species figure set and all plot source/audit tables.

```python
generate_mix_species_report(dsd: 'DrugSeqData', output_dir: 'str | Path', *, formats=('pdf', 'png'), expected_col: 'str | None' = 'compound', species_a_prefix: 'str' = 'ENSG', species_b_prefix: 'str' = 'ENSMUSG', species_a_name: 'str' = 'human', species_b_name: 'str' = 'mouse', expected_species_map: 'Mapping[str, str] | None' = None, assignment_threshold: 'float' = 0.5, min_species_umi: 'float' = 0, human_prefix: 'str | None' = None, mouse_prefix: 'str | None' = None, nrow_plate: 'int | None' = None, barnyard_metric_species: 'str | Sequence[str] | None' = 'base', palette: 'PaletteSpec' = None, expected_palette: 'PaletteSpec' = None, observed_palette: 'PaletteSpec' = None, plate_palette: 'PaletteSpec' = None, reference_line_color: 'str' = '#7F7F7F', guide_line_color: 'str' = '#CBD5E1', annotation_color: 'str' = '#D32F2F', edge_color: 'str' = 'white', grid_color: 'str' = '#DCE6F2', panel_border_color: 'str' = '#334155', point_color: 'str' = '#253746', label_color: 'str' = 'white', point_alpha: 'float' = 0.72, fill_alpha: 'float' = 0.95, show_observed_style: 'bool' = True) -> 'list[str]'
```

| Parameter | Type | Default | Meaning and use |
| --- | --- | --- | --- |
| `dsd` | DrugSeqData | required | Merged DrugSeqData input. This workflow does not depend on QC-report output. |
| `output_dir` | str \| Path | required | Directory receiving figures and source tables. |
| `formats` | — | ('pdf', 'png') | Iterable of figure formats; this HTML tutorial requests PNG only. |
| `expected_col` | str \| None | 'compound' | Optional metadata column containing the assigned base species, or labels that can be resolved with `expected_species_map`. |
| `species_a_prefix` | str | 'ENSG' | Gene-ID prefix selecting species A features. |
| `species_b_prefix` | str | 'ENSMUSG' | Gene-ID prefix selecting species B features. |
| `species_a_name` | str | 'human' | Display and assignment name for species A. |
| `species_b_name` | str | 'mouse' | Display and assignment name for species B. |
| `expected_species_map` | Mapping[str, str] \| None | None | Optional mapping from values in `expected_col` to species A, species B, `control`, or `unknown`. |
| `assignment_threshold` | float | 0.5 | Minimum within-pair fraction required for a species call. Use a value above 0.5 for mutually exclusive majority criteria; at a 0.5 tie, species A is selected first. This is a classification threshold, not a validated purity cutoff. |
| `min_species_umi` | float | 0 | Minimum species-A-plus-species-B UMI count; lower-count wells are called `low_umi`. |
| `human_prefix` | str \| None | None | Backward-compatible alias for `species_a_prefix`. |
| `mouse_prefix` | str \| None | None | Backward-compatible alias for `species_b_prefix`. |
| `nrow_plate` | int \| None | None | Maximum number of plate panels per row in the barnyard, purity, and concordance figures. |
| `barnyard_metric_species` | str \| Sequence[str] \| None | 'base' | Species selection forwarded to the barnyard `metric_species` option; defaults to assigned base species present in the data. |
| `palette` | PaletteSpec | None | Backward-compatible fallback palette for every categorical layer. Role-specific palettes below override it. |
| `expected_palette` | PaletteSpec | None | Colors for expected labels or cell types in the barnyard and purity plots. Accepts a Seaborn name, color sequence, category mapping, or `None`. |
| `observed_palette` | PaletteSpec | None | Fill colors for observed-species categories in the concordance bars. |
| `plate_palette` | PaletteSpec | None | Fill and point colors for plate groups in QC violins and species-fraction plots. |
| `reference_line_color` | str | '#7F7F7F' | Color of the barnyard balance diagonal and purity diagonal. |
| `guide_line_color` | str | '#CBD5E1' | Color of the 50% horizontal and vertical guides in the purity plot. |
| `annotation_color` | str | '#D32F2F' | Color of external barnyard summary statistics and QC median labels. |
| `edge_color` | str | 'white' | Point and stacked-bar edge color across the report. |
| `grid_color` | str | '#DCE6F2' | Color of dashed panel grid lines across all figures. |
| `panel_border_color` | str | '#334155' | Color of the complete border around every plot panel. |
| `point_color` | str | '#253746' | Color of individual points overlaid on the mix-species QC violins. |
| `label_color` | str | 'white' | Color of count labels placed inside concordance bars. |
| `point_alpha` | float | 0.72 | Opacity of barnyard/purity points and the fraction-plot point overlay; must be between zero and one. |
| `fill_alpha` | float | 0.95 | Opacity of stacked-bar and boxplot fills; must be between zero and one. |
| `show_observed_style` | bool | True | Whether barnyard and purity marker shapes also encode observed species. Set `False` for the cleaner reference-dashboard style. |

**Returns / side effects:** A list of figure paths plus CSV tables: per-well metrics, summaries, concordance counts, contamination-by-base data, audit candidates, and a figure manifest.

Set `expected_col='expected_species'` after `assign_well_species()` to make the intended species and off-target direction explicit and reproducible.


## Comparison, batch effects and reproducibility

Use these functions in the [single-plate](single-plate-comparison.md), [multi-plate](multi-plate-comparison.md), and [well-repeatability](well-to-well-repeatability.md) chapters. The [bulk comparison script](bulk-rnaseq-comparison.md#run) is a separate tutorial helper, not a package API.

### `normalize_counts`

```python
ds.normalize_counts(dsd, method='limma_voom', log_transform=True, prior_count=1.0, target_sum=1000000.0, inplace=True)
```

Replaces `adata.X` with normalized expression and preserves `layers['counts']`. The `limma_voom` option produces TMM-normalized log2-CPM; design-specific precision weights are estimated during DE. `prior_count` applies to the CPM/TMM transformations; limma-voom uses its fixed count offset. Returns `None` in place or a copied object with `inplace=False`.

### `plot_comparison_heatmap`

```python
ds.plot_comparison_heatmap(dsd, contrasts=None, n_top_genes=60, value_col='logFC', fdr_threshold=0.05, lfc_threshold=0.5, significant_only=True, cluster_rows=True, cluster_cols=True, cmap='RdBu_r', figsize=None)
```

Returns a Matplotlib figure containing gene-level effects across contrasts. The gene list is the union of the top absolute-effect genes per selected contrast, significant only by default. This differs from a sample-by-gene expression heatmap and from a well-by-well correlation heatmap. Save the figure with `fig.savefig(...)`.


These signatures describe the installed DRUGseqPy 0.1.0 implementation used to verify this tutorial. Positional parameters and defaults are preserved; type annotations are omitted for readability.

### `audit_control_outliers`

```python
ds.audit_control_outliers(dsd, control_label='DMSO', control_col='compound', score_col='outlier_score', threshold=5.0, columns=None)
```

Returns a control-only pandas DataFrame. Does not modify or remove wells.

The returned `excluded` flag means the score exceeds the supplied threshold; apply exclusions explicitly only after review. The default threshold is illustrative, not a universal QC cutoff.

### `plot_control_outlier_audit`

```python
ds.plot_control_outlier_audit(audit, x='total_umi', y='n_genes_det', label_col='well_id', plate_col='plate_id', figsize=(8.3, 5.8))
```

Returns a Matplotlib Figure; does not change the data.

Its legend says “Excluded control” for flagged rows even though the audit has not removed them. Clarify this in captions or rename the display labels.

### `correct_batch`

```python
ds.correct_batch(dsd, batch_col='plate_id', method='limma_voom', covariate_cols=None, prior_count=1.0, inplace=True)
```

With `inplace=True`, returns None and replaces X with corrected log2-CPM; with False, returns an independent DrugSeqData. Saves batch_uncorrected and batch_corrected layers and batch settings in uns.

Raw counts remain unchanged. Preserve biological covariates and check the design for confounding. For a fair before/after comparison, take the before view from batch_uncorrected. The edgeR option here is a Python TMM expression transformation; differential-expression edgeR separately needs R/edgeR and rpy2.

### `compute_batch_effect_metrics`

```python
ds.compute_batch_effect_metrics(before, after, batch_col='plate_id', group_col='compound', reduction='X_pca', n_dims=10, random_state=42)
```

Returns a pandas DataFrame with metric, before and after columns. Does not mutate the inputs.

Both objects must have the requested embedding and exactly the same samples in the same order. Plate mixing alone does not prove that biology has been preserved.

### `run_pca`

```python
ds.run_pca(dsd, n_pcs=50, n_variable_genes=2000, use_highly_variable=True, random_state=42, inplace=True)
```

With inplace=True, returns None; otherwise returns a copy. Writes X_pca to obsm, PCs to varm, PCA statistics to uns, and highly_variable to var.

Normalize X before use. Use use_highly_variable=False for the TMM/log2-CPM comparison shown here. The number of components is capped by sample and feature counts.

### `run_umap`

```python
ds.run_umap(dsd, dims=20, n_neighbors=15, min_dist=0.3, metric='euclidean', use_rep='X_pca', random_state=42, inplace=True)
```

With inplace=True, returns None; otherwise returns a copy. Writes X_umap and the Scanpy neighbor graph.

Requires a precomputed representation, normally X_pca. A fixed random_state helps reproducibility; the embedding is exploratory.

### `run_tsne`

```python
ds.run_tsne(dsd, dims=20, perplexity=30.0, use_rep='X_pca', random_state=42, inplace=True)
```

With inplace=True, returns None; otherwise returns a copy. Writes X_tsne and parameters in uns.

Requires a precomputed representation and at least three samples. Perplexity is capped at n_samples minus one.

### `compute_within_plate_repeatability`

```python
ds.compute_within_plate_repeatability(dsd, plate_col='plate_id', group_col='compound', min_replicates=2, use_norm=False, method='pearson')
```

Returns a pandas DataFrame of within-plate, within-treatment replicate-pair correlations. Does not mutate the object.

use_norm=False uses log10(raw UMI + 1); True uses X as supplied. Correlations between replicate wells measure repeatability, not treatment-effect reproducibility.

### `summarise_within_plate_repeatability`

```python
ds.summarise_within_plate_repeatability(repeatability, plate_col='plate_id', group_col='compound', value_col='correlation')
```

Returns a pandas DataFrame summarizing the pair table.

Supply the pair table from compute_within_plate_repeatability. Pairs sharing a well are statistically dependent.

### `compute_cross_plate_reproducibility`

```python
ds.compute_cross_plate_reproducibility(dsd, plate_col='plate_id', group_col='compound', reference='DMSO', exclude_levels=('water',), method='limma_voom', min_replicates=2, fdr_threshold=0.05, lfc_threshold=0.5)
```

Returns a dictionary containing summary, plate_results, gene_results and config. Does not store these results in the input object.

Fits DE separately on each plate. Every requested treatment needs adequate replication and controls on each plate. The normalized-expression fields use X as supplied, so normalize first. Gene-level correlations can be undefined for constant profiles.

### `run_comparison`

```python
ds.run_comparison(dsd, case, control, group_col='compound', method='limma_voom', batch_cols=None, within_plate=True, min_replicates=2, fdr_threshold=0.05, lfc_threshold=0.5, contrast_name=None)
```

Returns a DE pandas DataFrame for one contrast; uses the raw counts layer.

Positive logFC means higher expression in case than control. Include batch_cols for a plate-adjusted model. within_plate restricts controls to plates represented among the case samples; it does not itself fit separate per-plate models.

### `compute_multi_de`

```python
ds.compute_multi_de(dsd, compounds=None, reference='DMSO', compound_col='compound', method='limma_voom', within_plate=True, min_replicates=2, n_jobs=1, fdr_threshold=0.05, lfc_threshold=0.5, inplace=True, *, methods=None, comparisons=None, group_col=None, batch_cols=None, exclude_levels=None)
```

With inplace=True, returns None; otherwise returns an independent DrugSeqData. Rebuilds comparison_results and comparison_failures in uns; default compound-versus-reference runs also populate de_results.

Check uns["comparison_failures"] because failed contrasts can be recorded as warnings. Raw counts supply DE. Multiple requested methods use keys such as contrast::method.

### `summarise_comparison`

```python
ds.summarise_comparison(dsd, fdr_threshold=0.05, lfc_threshold=0.5)
```

Returns a pandas DataFrame of contrast-level tested/significant gene counts and cutoffs.

Check for an empty table before choosing a priority contrast. Selecting the strongest example is a presentation choice, not independent validation.

### `run_enrichment`

```python
ds.run_enrichment(dsd, contrasts=None, libraries=('hallmark', 'go_bp', 'reactome', 'kegg'), gene_set_paths=None, mode=('gsea', 'ora'), rank_by='stat', fdr_threshold=0.05, lfc_threshold=0.5, n_perm=1000, min_size=15, max_size=500, species='Human', background='all_tested', seed=42, cache_dir=None, allow_network=False, inplace=True)
```

With inplace=True, returns None; otherwise returns an independent DrugSeqData. Stores per-contrast libraries in enrichment_results and the flattened summary in enrichment_summary.

Provide matching local GMT files for offline runs; allow_network=False is the default. Match organism and gene identifiers. Directional ORA uses significant up/down genes and an explicit tested-gene background. GSEA uses ranked genes.

### `enrichment_summary`

```python
ds.enrichment_summary(dsd)
```

Returns the flattened enrichment pandas DataFrame.

An empty table can mean there are no eligible gene sets, no significant genes in a direction, or no stored enrichment results; inspect run diagnostics.

### `run_comparison_workflow`

```python
ds.run_comparison_workflow(dsd, output_dir=None, *, comparisons=None, group_col='compound', reference='DMSO', exclude_levels=None, methods=('limma_voom',), batch_cols=None, within_plate=True, min_replicates=2, fdr_threshold=0.05, lfc_threshold=0.5, n_jobs=1, normalize=None, normalize_kwargs=None, filter_samples_kwargs=None, filter_genes_kwargs=None, enrichment=True, libraries=('hallmark', 'go_bp', 'reactome', 'kegg'), gene_set_paths=None, enrichment_mode=('gsea', 'ora'), n_perm=1000, min_size=15, max_size=500, species='Human', background='all_tested', seed=42, cache_dir=None, allow_network=False, figures=('volcano', 'heatmap', 'bar', 'box', 'violin', 'enrichment'), n_top_genes=60, boxplot_genes=None, formats=('png', 'pdf'), inplace=True)
```

Returns a dictionary with summary, manifest, figures, files and config, plus dsd when inplace=False. With an output directory, exports result tables, figures and configuration.

Enrichment defaults to True: set enrichment=False for a minimal offline run. Inspect manifest statuses failed/partial and ensure summary is nonempty. It catches some stage failures instead of raising. Filtering is opt-in; normalize changes X while DE still uses raw counts.

### Plotting helpers

These functions return a Matplotlib Figure and do not save it automatically. Use `fig.savefig(...)`, create its parent directory first, and close the figure when finished. They consume the objects or result tables described above.

```python
ds.plot_pca(dsd, color_by='compound', **kwargs)
```

```python
ds.plot_umap(dsd, color_by='compound', **kwargs)
```

```python
ds.plot_tsne(dsd, color_by='compound', **kwargs)
```

```python
ds.plot_mds(dsd, group_by='sample_type', split_by=None, label_by='compound', n_top_genes=500, use_norm=True, figsize=None)
```

```python
ds.plot_within_plate_repeatability(repeatability, plate_col='plate_id', group_col='compound', value_col='correlation', figsize=(14, 6.2))
```

```python
ds.plot_within_plate_replicate_scatter(dsd, plate, compound, method='pearson', plate_col='plate_id', group_col='compound', label_col='well_id', use_norm=False, max_pairs=4, figsize=None)
```

```python
ds.plot_within_plate_correlation_heatmap(dsd, plate, compound, method='pearson', plate_col='plate_id', group_col='compound', label_col='well_id', use_norm=False, figsize=(7.2, 6.2))
```

```python
ds.plot_single_plate_de_volcano(reproducibility, plate, compound, fdr_threshold=0.05, lfc_threshold=0.5, n_label=5, figsize=(7, 5.8))
```

```python
ds.plot_single_plate_de_layout(dsd, reproducibility, plate, compound, plate_col='plate_id', n_top_genes=30, nrow=8, ncol=12)
```

```python
ds.plot_single_plate_de_heatmap(dsd, reproducibility, plate, compound, plate_col='plate_id', group_col='compound', n_top_genes=30, use_norm=True, figsize=(11, 7))
```

```python
ds.plot_cross_plate_reproducibility(reproducibility, correlation_threshold=0.5, correlation_method='pearson', figsize=(13, 6.2))
```

```python
ds.plot_plate_logfc_concordance(reproducibility, compounds=None, max_panels=4, value_type='log2fc', correlation_method='pearson', figsize=(11.5, 10.5))
```

```python
ds.plot_deg_counts(dsd, contrasts=None, fdr_threshold=0.05, lfc_threshold=0.5, figsize=None)
```

```python
ds.plot_comparison_volcano(dsd, contrast=None, method=None, lfc_threshold=0.5, fdr_threshold=0.05, n_label=15, figsize=(6, 5))
```

```python
ds.plot_enrichment_dotplot(dsd, contrast=None, libraries=None, direction='up', mode='ora', n_terms=20, figsize=None)
```

For `plot_plate_logfc_concordance`, `value_type="normalized"` uses the X values present when reproducibility was computed. Use `value_type="log2fc"` to compare treatment effects. Plotting cutoffs should match the DE run.

---

[← Previous chapter](bulk-rnaseq-comparison.md) · [Next chapter →](troubleshooting.md)
