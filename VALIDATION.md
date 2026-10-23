# Validation

The tutorial package was checked on 23 October 2026. The analysis environment used Python 3.10 and a local drugseqpy development build reporting version 0.1.0. The documentation was built with Python 3.12. This does not establish that a particular published package revision provides the same implementation.

## Analysis checks

- The original bundled synthetic run completed QC, mixed-species, and multi-plate comparison. The expanded `scripts/run_tutorial.py all` also completed the new single-plate workflow, including a run with reviewed exclusions and an offline synthetic GMT.
- The original three-workflow run produced 20 figures and 30 CSV/TSV tables. Both treatment contrasts succeeded with no failed or partial comparison-manifest entries.
- The reviewed-exclusion option completed a separate comparison run: the requested well was removed, its reason preserved, control counts updated, and both contrasts succeeded.
- All four bundled metadata tables pass `validate_metadata`. The seed, data classification, and checksums are in `data/practice-manifest.json`.
- All 23 non-STAR QC blocks, 13 mixed-species blocks, and 48 non-enrichment comparison blocks completed sequentially against the bundled inputs. This included both PCA, UMAP, and t-SNE views, reproducibility, and differential-expression plots.
- Input-guide examples were exercised with shuffled metadata, terminal barcode suffixes, a metadata-only well, duplicate-ID prevention, and physical-coordinate reconciliation.
- The 37 comparison figures were regenerated through the worked pipeline with revised labels. That run exercised PCA, UMAP, t-SNE, differential expression, reproducibility, and enrichment.
- The 31 QC/species figures had their sample-label text revised from vector figures, with plotted data preserved. Figure metadata were stripped, extracted labels checked, and all figure groups visually reviewed.

- The single-plate practice route retained 36 biological wells and 800 genes and fitted both treatment contrasts. Its saved configuration uses no plate covariate and `within_plate=False`. All 15 non-enrichment chapter code blocks ran sequentially; the chapter's stronger gene filter retained 795 genes.
- The bulk practice helper compared 16 plate/reference profiles and 368 well/reference pairs on 400 shared genes per species. Checks covered version suffixes, missing measurements, duplicate stable IDs, negative values, constant vectors, incompatible conditions, reviewed exclusions, and raw-count preservation. Its independent synthetic generator reproduced its manifest checksums.
- The well-repeatability helper computed all 264 matched-condition pairs on 48 synthetic wells. Its coefficients match independent SciPy calculations; pair/heatmap gene masks, undefined constant profiles, explicit pair selection, retained flags, sorting alignment, and display-only jitter were checked. Five synthetic illustrative figures were added.
- Nine single-plate figures were added. The four QC panels were regenerated with identical QC values and neutral labels; five DE/enrichment panels retain their original pixels and numerical results. A rerun on the current development build changes some DE estimates, so the historical snapshots are not represented as exact current-version outputs.
- Fourteen bulk comparison figures were anonymized by editing vector text, preserving plotted values and geometry. The six historical evidence tables retain all 160 body cells after label substitution. Reference accession labels were removed from the public figures and prose.

The default bundled workflow and the detailed chapter examples passed. Enrichment with biological gene-set collections and STARsolo read-level QC require separate optional inputs. The small synthetic GMT checks validate the enrichment interface, not the historical pathway findings.

## Documentation checks

All eleven chapters passed browser checks at 1440 px and 390 px: no horizontal page overflow, chapter navigation present, and mobile navigation opens and closes. Desktop and mobile screenshots were visually reviewed. Every figure has a working link to its full-size asset. All 96 PNG files decode successfully; the static checker verifies their references.


The expanded static validation passes 32 documents, 1,073 local links/assets, 115 Python block syntax checks, all 96 figure references, and both synthetic-data checksum manifests. Run `python scripts/validate_site.py` after rebuilding. It checks local page/asset links, fragment targets, duplicate HTML IDs, runnable Python-block syntax, image descriptions, old sample identifiers, private filesystem references, and synthetic-data checksums. API signature displays are not executable examples.

The eleven Markdown chapters are the editable source for eleven generated HTML chapters. Six old HTML names redirect to the corresponding revised chapters. Figures and browser assets are local, and all internal site links are relative so the site works below a GitHub Pages repository path.

## Reproduction boundaries

The worked figures and tables are fixed analysis examples with Sample 1–4 labels. The bundled practice data were generated independently and do not reproduce those numerical values. The practice runner deliberately uses a smaller gene-filtering rule, explicit species review settings, and no automatic DMSO exclusions. Bulk practice outputs are descriptive and do not recreate the historical bootstraps or pseudogene audit. Chapter text explains these differences.

Some report functions warn when optional housekeeping genes are not present in the synthetic feature list; the warning does not invalidate the core QC output. Matching the displayed package version alone is not enough to reproduce an analysis. Record the exact compatible source revision, environment, inputs, parameters, exclusions, and gene-set versions for real work.
