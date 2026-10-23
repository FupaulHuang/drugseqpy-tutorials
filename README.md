# drugseqpy tutorials

A connected beginner tutorial for count-matrix preparation, quality control, mixed-species analysis, well-to-well repeatability, single-plate and multi-plate treatment comparisons, and bulk RNA-seq validation.

**[Visit the tutorial website](https://fupaulhuang.github.io/drugseqpy-tutorials/)** · **[Read on GitHub](docs/index.md)** · **[Get started](docs/getting-started.md)**

The website presents the connected chapters as readable tutorial pages. The same content is available as Markdown in this repository. To browse offline, open `docs/index.html` locally.

The chapters retain anonymized worked results with Sample 1–4 labels. Four separate, clearly labeled **synthetic practice datasets** and a portable runner make the workflow runnable without the worked inputs.

## Read in order

1. [Getting started](docs/getting-started.md)
2. [Prepare inputs](docs/input-preparation.md)
3. [Quality control](docs/qc-workflow.md)
4. [Mixed-species analysis](docs/mix-species-workflow.md)
5. [Well-to-well repeatability](docs/well-to-well-repeatability.md)
6. [Single-plate comparison](docs/single-plate-comparison.md)
7. [Multi-plate comparison](docs/multi-plate-comparison.md)
8. [Bulk RNA-seq comparison](docs/bulk-rnaseq-comparison.md)
9. [API reference](docs/api-reference.md)
10. [Troubleshooting](docs/troubleshooting.md)

Sample 1 is the single-plate treatment example; Sample 1 and Sample 2 together demonstrate multi-plate comparison. Sample 3 and Sample 4 are mixed-species plate datasets. Each observation inside a dataset is one well.

## Run the practice workflow

Install a compatible drugseqpy build following [Getting started](docs/getting-started.md), then run from this repository's root:

```bash
python scripts/check_environment.py
python scripts/run_tutorial.py all --data-dir data --output-dir results
```

The command above runs QC, species review, single-plate comparison, and multi-plate comparison. The [repeatability](docs/well-to-well-repeatability.md) and [bulk RNA-seq](docs/bulk-rnaseq-comparison.md) chapters provide separate runnable scripts.

The examples were checked against a local development build reporting version 0.1.0. Public revisions may expose different APIs; the environment check verifies required functions and parameters. A matching public release is not assumed. The practice runner uses PCA by default; optional UMAP/t-SNE and enrichment are described in `--help`.

## Preview the website

```bash
python -m http.server 8000 --directory docs
```

Visit `http://localhost:8000`. All styles, scripts, and figures are local; there are no external fonts, analytics, or embedded data downloads.

## GitHub Pages deployment

This repository publishes the `docs/` folder on the `main` branch at [fupaulhuang.github.io/drugseqpy-tutorials](https://fupaulhuang.github.io/drugseqpy-tutorials/). Pushing updated HTML and assets to that branch triggers a site update. `docs/.nojekyll` disables Jekyll processing.

To publish a fork, open **Settings → Pages**, choose **Deploy from a branch**, select your branch and **/docs**, and save. See the [GitHub Pages publishing guide](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site). The six previous HTML filenames remain as redirect pages, so existing chapter links still lead into the connected site.

## Maintain the tutorials

The Markdown chapters in `docs/` are the editable source. Rebuild the matching HTML after changes:

```bash
python -m pip install -r requirements-docs.txt
python scripts/build_site.py
python scripts/validate_site.py
```

The build does not need any experimental inputs. The included figures are curated static assets. To create a fresh synthetic input copy, run `python scripts/make_practice_data.py --output practice-copy` in the analysis environment. It refuses to replace existing sample folders.

| Location | Contents |
| --- | --- |
| `docs/*.md` | Complete Markdown chapters |
| `docs/*.html` | Ready-to-publish website |
| `docs/assets/` | Shared styling, browser behavior, and anonymized figures |
| `data/` | Synthetic practice matrices, metadata, and checksums |
| `scripts/` | Practice-data generator, workflow runner, environment check, site build, and validation |
| `VALIDATION.md` | Scope and results of the release checks |
