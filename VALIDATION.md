# Validation

The tutorial package was checked on 23 October 2026. The analysis environment used Python 3.10 and a local drugseqpy development build reporting version 0.1.0. The documentation was built with Python 3.12. This does not establish that a particular published package revision provides the same implementation.

## Analysis checks

- The bundled synthetic inputs completed the QC, mixed-species, and comparison workflows with `scripts/run_tutorial.py all`.
- The run produced 20 figures and 30 CSV/TSV tables. Both treatment contrasts succeeded with no failed or partial comparison-manifest entries.
- The reviewed-exclusion option completed a separate comparison run: the requested well was removed, its reason preserved, control counts updated, and both contrasts succeeded.
- All four bundled metadata tables pass `validate_metadata`. The seed, data classification, and checksums are in `data/practice-manifest.json`.
- All 23 non-STAR QC blocks, 13 mixed-species blocks, and 48 non-enrichment comparison blocks completed sequentially against the bundled inputs. This included both PCA, UMAP, and t-SNE views, reproducibility, and differential-expression plots.
- Input-guide examples were exercised with shuffled metadata, terminal barcode suffixes, a metadata-only well, duplicate-ID prevention, and physical-coordinate reconciliation.
- The 37 comparison figures were regenerated through the worked pipeline with revised labels. That run exercised PCA, UMAP, t-SNE, differential expression, reproducibility, and enrichment.
- The 31 QC/species figures had their sample-label text revised from vector figures, with plotted data preserved. Figure metadata were stripped, extracted labels checked, and all figure groups visually reviewed.

The default bundled workflow and the detailed chapter examples passed. Enrichment and STARsolo read-level QC require separate optional inputs and are not claimed as bundled-data tests.

## Documentation checks

All eight chapters passed browser checks at 1440 px and 390 px: no horizontal page overflow, chapter navigation present, and mobile navigation opens and closes. Desktop and mobile screenshots were visually reviewed. All 68 PNG files decode successfully; the static checker verifies their references.


Run `python scripts/validate_site.py` after rebuilding. It checks local page/asset links, fragment targets, duplicate HTML IDs, runnable Python-block syntax, image descriptions, old sample identifiers, private filesystem references, and synthetic-data checksums. API signature displays are not executable examples.

The eight Markdown chapters are the editable source for eight generated HTML chapters. Four old HTML names redirect to the corresponding revised chapters. Figures and browser assets are local, and all internal site links are relative so the site works below a GitHub Pages repository path.

## Reproduction boundaries

The worked figures and tables are fixed analysis examples with Sample 1–4 labels. The bundled practice data were generated independently and do not reproduce those numerical values. The practice runner deliberately uses a smaller gene-filtering rule, explicit species review settings, and no automatic DMSO exclusions. Chapter text explains these differences.

Some report functions warn when optional housekeeping genes are not present in the synthetic feature list; the warning does not invalidate the core QC output. Matching the displayed package version alone is not enough to reproduce an analysis. Record the exact compatible source revision, environment, inputs, parameters, exclusions, and gene-set versions for real work.
