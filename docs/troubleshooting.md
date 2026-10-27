# Troubleshooting

[Learning path](index.md) · [Input checks](input-preparation.md) · [API reference](api-reference.md)

## A tutorial function is missing

Run `python scripts/check_environment.py` in the same environment used for analysis. The examples need a recent development API, including comparison and report functions. The package version can remain `0.1.0` across revisions. Install a compatible checkout with `python -m pip install /path/to/drugseqpy` if the available published revision does not provide the functions. The tutorials were validated against a local development build; a matching public release is not assumed.

## Matrix or barcode loading fails

Pass the `matrix/` directory containing all three **gzipped** files directly. Use sparse coordinate Matrix Market format, with features in matrix-row order and barcodes in matrix-column order. Metadata can be in another row order because import joins by identifier. Specify `metadata_key="barcode"` when a metadata table has both barcode and sample ID columns. Compare the identifiers after terminal `-1` removal and inspect duplicates before import.

See [Prepare inputs](input-preparation.md) for executable checks. `min_overlap=1.0` ensures every matrix barcode matches metadata; it does not imply every planned metadata well appears in a filtered matrix.

## A plate map looks transposed

Parseable `well_id` values take precedence over numeric coordinates. Well-dT numbering runs down the plate column. Confirm the physical scheme, keep original well labels in a separate column, and use verified A1-style positions on a plotting copy. Check more than one corner and compare the generated layout against the design. See the [coordinate recipe](input-preparation.md#physical-coordinates).

## QC reports a missing dose

A species-control experiment may have no meaningful treatment dose. Decide whether missing dose is intentional and record it; do not invent a numeric dose to make an audit pass. A treatment comparison needs unambiguous dose and time labels. Compare matched conditions and use a joint condition label when compound alone would pool different doses or times.

## Species calls are unexpected

Keep stable species-prefixed IDs active in `var_names`. Confirm both species prefixes match the reference and that the expected-species mapping covers every design label. A human/mouse gene symbol alone does not establish species identity. Species fractions use only the two configured species in their denominator; inspect unassigned features separately. Low counts can make fractions unstable, and a majority call at 0.5 is not a stringent purity requirement. Pass identical thresholds to tables and figures.

## There are no STARsolo read-level figures

The synthetic practice bundles contain count matrices, not read alignments or read statistics. Skip the optional STARsolo section unless you supply the matching `Summary.csv`, `CellReads.stats`, and barcode structure. Raw reads, UMI counts, and normalized counts answer different questions. A read-depth sensitivity analysis requires actual downsampling and reprocessing; library-size normalization does not equalize sequencing reads.

## I have only one plate

Use the [single-plate comparison](single-plate-comparison.md) and run `python scripts/run_tutorial.py single-plate`. It compares treatments against controls on Sample 1 with `batch_cols=None` and `within_plate=False`. One plate cannot establish between-plate effects or cross-plate treatment-effect concordance.

## Batch correction fails or removes treatment separation

Review treatment × plate counts after every exclusion. Cases on one plate and controls on another are confounded. Preserve estimable biological terms in `covariate_cols`, and use the same normalization for before/after views. Inspect treatment structure as well as plate mixing. Use raw counts with a plate term for differential expression, even when exploratory plots use corrected expression.

## Comparison returns no results

Check exact DMSO spelling, sample counts per condition and plate, the reference group, feature filters, and the saved comparison manifest. A successful call can record unsupported contrasts in its manifest. An empty significant-gene list is a valid outcome; it differs from a failed or skipped model. Do not lower thresholds solely to obtain significant genes.

## Enrichment is empty or GMT files are missing

The first-run practice workflow skips enrichment. For an optional enrichment run, provide local GMT paths whose species and identifiers match the tested genes. Check gene-set overlap and the tested-gene background. An empty upregulated or downregulated result can be valid. Keep the library version and analysis parameters with your results. No gene-set download is performed by default.

## Which term should I use for repeated measurements?

Use **within-plate repeatability** for matched wells, **between-plate consistency within one batch** for the measured Drug-seq case, and **cross-plate treatment-effect concordance** for correlations of separately fitted log2 fold changes. Experimental reproducibility requires explicitly changed conditions; merely fitting several plates does not demonstrate it. API names containing `reproducibility` remain valid function names. [See the terminology guide](input-preparation.md#repeatability-and-concordance).

## Numbers differ from the figures

The detailed chapters show fixed anonymized worked results. The bundled practice counts are synthetic and intentionally produce different numbers. Your own data will differ too. Compare table schemas, analysis decisions, and interpretation. For numerical reproducibility, record the exact package revision, dependencies, filters, exclusions, seeds, and input checksums; the package version alone is insufficient.

## A scatter correlation differs from its heatmap cell

Check which genes enter each calculation. The [well-repeatability workflow](well-to-well-repeatability.md) uses the plate-wide detected-gene set for the full heatmap, but the union of genes detected in the two selected wells for each scatter. Zero-only genes in that pair are excluded from the scatter correlation, so the coefficients can differ. Match the transformation and gene universe before comparing values.

## Bulk agreement changes across reference studies

Review cell identity, condition, library type, expression units, gene intersections, and whether each reference was declared primary or sensitivity evidence. Pair Sample 3 and Sample 4 against the same reference and summarize within studies before combining them. The [bulk comparison chapter](bulk-rnaseq-comparison.md) explains why correlations and conditional bootstrap intervals do not remove plate–method confounding.

## Read the website locally or publish it

Open `docs/index.html` directly, or run `python -m http.server 8000 --directory docs` from the repository root and visit `http://localhost:8000`. No network connection is needed to read the site. Code-copy buttons need browser clipboard support; if unavailable, they select the code for manual copying.

For GitHub Pages, follow the publishing instructions in the repository README. Both Markdown and HTML are included; edit the Markdown and run `python scripts/build_site.py` to keep them synchronized.

[← API reference](api-reference.md) · [Return to the learning path](index.md)
