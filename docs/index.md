# Learn drugseqpy

Build a complete analysis from count matrices and plate metadata: inspect well quality and species identity, assess well-to-well repeatability, compare treatments within and across plates, and assess agreement with bulk RNA-seq.

These tutorials are written for new users. Each chapter explains the inputs, shows the Python calls, and connects the outputs to the next decision. You can read the same chapters as a website or directly on GitHub.

## Follow the learning path

| Step | Chapter | What you will be able to do |
| --- | --- | --- |
| 1 | [Getting started](getting-started.md) | Install the package and run a small practice workflow. |
| 2 | [Prepare inputs](input-preparation.md) | Align barcodes, validate counts and metadata, and check physical well positions. |
| 3 | [Quality control](qc-workflow.md) | Review library metrics, plate maps, marker expression, and optional read-level QC. |
| 4 | [Mixed-species analysis](mix-species-workflow.md) | Assign intended species and inspect purity, discordance, and reciprocal off-target signal. |
| 5 | [Well-to-well repeatability](well-to-well-repeatability.md) | Read well-correlation heatmaps and inspect treatment-matched replicate scatter plots. |
| 6 | [Single-plate comparison](single-plate-comparison.md) | Audit controls and estimate treatment effects within one plate. |
| 7 | [Multi-plate comparison](multi-plate-comparison.md) | Audit controls, diagnose batch effects, assess treatment-effect concordance, and fit plate-adjusted treatment effects. |
| 8 | [Drug-seq worked case](drugseq-case-study.md) | Analyze a selected real Drug-seq experiment from one plate to three matched plates. |
| 9 | [Bulk RNA-seq comparison](bulk-rnaseq-comparison.md) | Align features, compare matched reference profiles, and respect study and plate structure. |
| 10 | [API reference](api-reference.md) | Find function signatures, parameters, returns, and related workflow steps. |

Keep [Troubleshooting](troubleshooting.md) nearby when adapting the examples.

## Understand the four sample labels

| Dataset label | Role in the tutorial | Folder |
| --- | --- | --- |
| Sample 1 | Single-plate example and first plate in the multi-plate comparison | `data/sample_1/` |
| Sample 2 | Second plate with matching treatments and controls | `data/sample_2/` |
| Sample 3 | First mixed-species plate for QC and species review | `data/sample_3/` |
| Sample 4 | Second mixed-species plate for QC and species review | `data/sample_4/` |

The [Drug-seq worked case](drugseq-case-study.md) uses a **separate cohort**, with case-specific Sample 1–3 labels and its own input directory. Those labels do not refer to the datasets in this table.

Here **Sample 1–4 are labels for plate datasets**. Within a dataset, each observation is one well, and `sample_id` identifies that well. Sample 3 and Sample 4 form a different example from Sample 1 and Sample 2; do not pool all four for a treatment comparison.

The detailed chapters retain anonymized worked figures and result tables. The bundled **synthetic practice data** are independently generated and contain no measurements from the worked example. They teach the same input schema and analysis calls, but their numerical results will differ. The practice data include low-depth and mixed wells to make the review steps useful.

## Choose a route

Start with [Getting started](getting-started.md) to run the included practice data. If you already have count matrices, continue to [Prepare inputs](input-preparation.md), then run the same workflow on a separate data directory. Mixed-species analysis is optional for a single-species experiment. For treatment analysis, start with [Single-plate comparison](single-plate-comparison.md), then extend to [Multi-plate comparison](multi-plate-comparison.md) when matched treatments and controls are available on another plate. The [Drug-seq worked case](drugseq-case-study.md) applies the workflow to three selected real plates with a compact matched treatment panel. The [repeatability chapter](well-to-well-repeatability.md) helps review replicate consistency before treatment modeling. The [bulk comparison](bulk-rnaseq-comparison.md) extends Sample 3/4 QC to independently prepared references. Enrichment requires suitable local gene-set files; STARsolo read-level QC requires separate read statistics.

[Start the tutorial →](getting-started.md)
