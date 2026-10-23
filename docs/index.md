# Learn drugseqpy

Build a complete analysis from count matrices and plate metadata: inspect well quality, check mixed-species signal, compare treatments across plates, and interpret the saved results.

These tutorials are written for new users. Each chapter explains the inputs, shows the Python calls, and connects the outputs to the next decision. You can read the same chapters as a website or directly on GitHub.

## Follow the learning path

| Step | Chapter | What you will be able to do |
| --- | --- | --- |
| 1 | [Getting started](getting-started.md) | Install the package and run a small practice workflow. |
| 2 | [Prepare inputs](input-preparation.md) | Align barcodes, validate counts and metadata, and check physical well positions. |
| 3 | [Quality control](qc-workflow.md) | Review library metrics, plate maps, marker expression, and optional read-level QC. |
| 4 | [Mixed-species analysis](mix-species-workflow.md) | Assign intended species and inspect purity, discordance, and reciprocal off-target signal. |
| 5 | [Multi-plate comparison](multi-plate-comparison.md) | Audit controls, diagnose batch effects, measure reproducibility, and fit plate-adjusted treatment effects. |
| 6 | [API reference](api-reference.md) | Find function signatures, parameters, returns, and related workflow steps. |

Keep [Troubleshooting](troubleshooting.md) nearby when adapting the examples.

## Understand the four sample labels

| Dataset label | Role in the tutorial | Folder |
| --- | --- | --- |
| Sample 1 | First plate in the treatment comparison | `data/sample_1/` |
| Sample 2 | Second plate with matching treatments and controls | `data/sample_2/` |
| Sample 3 | First mixed-species plate for QC and species review | `data/sample_3/` |
| Sample 4 | Second mixed-species plate for QC and species review | `data/sample_4/` |

Here **Sample 1–4 are labels for plate datasets**. Within a dataset, each observation is one well, and `sample_id` identifies that well. Sample 3 and Sample 4 form a different example from Sample 1 and Sample 2; do not pool all four for a treatment comparison.

The detailed chapters retain anonymized worked figures and result tables. The bundled **synthetic practice data** are independently generated and contain no measurements from the worked example. They teach the same input schema and analysis calls, but their numerical results will differ. The practice data include low-depth and mixed wells to make the review steps useful.

## Choose a route

Start with [Getting started](getting-started.md) to run the included practice data. If you already have count matrices, continue to [Prepare inputs](input-preparation.md), then run the same workflow on a separate data directory. Mixed-species analysis is optional for a single-species experiment. Enrichment requires suitable local gene-set files; STARsolo read-level QC requires separate read statistics.

[Start the tutorial →](getting-started.md)
