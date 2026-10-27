# Validation

The tutorial package was checked on 23 October 2026. The analysis environment used Python 3.10 and a local drugseqpy development build reporting version 0.1.0. The documentation was built with Python 3.12. This does not establish that a particular published package revision provides the same implementation.

## Analysis checks

- The original bundled synthetic run completed QC, mixed-species, and multi-plate comparison. The expanded `scripts/run_tutorial.py all` also completed the new single-plate workflow, including a run with reviewed exclusions and an offline synthetic GMT.
- The original three-workflow run produced 20 figures and 30 CSV/TSV tables. Both treatment contrasts succeeded with no failed or partial comparison-manifest entries.
- The reviewed-exclusion option completed a separate comparison run: the requested well was removed, its reason preserved, control counts updated, and both contrasts succeeded.
- All four bundled metadata tables pass `validate_metadata`. The seed, data classification, and checksums are in `data/practice-manifest.json`.
- All 23 non-STAR QC blocks, 13 mixed-species blocks, and 48 non-enrichment comparison blocks completed sequentially against the bundled inputs. This included both PCA, UMAP, and t-SNE views, cross-plate treatment-effect concordance, and differential-expression plots.
- Input-guide examples were exercised with shuffled metadata, terminal barcode suffixes, a metadata-only well, duplicate-ID prevention, and physical-coordinate reconciliation.
- The 37 comparison figures were regenerated through the worked pipeline with revised labels. That run exercised PCA, UMAP, t-SNE, differential expression, cross-plate treatment-effect concordance, and enrichment.
- The 31 QC/species figures had their sample-label text revised from vector figures, with plotted data preserved. Figure metadata were stripped, extracted labels checked, and all figure groups visually reviewed.

- The single-plate practice route retained 36 biological wells and 800 genes and fitted both treatment contrasts. Its saved configuration uses no plate covariate and `within_plate=False`. All 15 non-enrichment chapter code blocks ran sequentially; the chapter's stronger gene filter retained 795 genes.
- The bulk practice helper compared 16 plate/reference profiles and 368 well/reference pairs on 400 shared genes per species. Checks covered version suffixes, missing measurements, duplicate stable IDs, negative values, constant vectors, incompatible conditions, reviewed exclusions, and raw-count preservation. Its independent synthetic generator reproduced its manifest checksums.
- The well-repeatability helper computed all 264 matched-condition pairs on 48 synthetic wells. Its coefficients match independent SciPy calculations; pair/heatmap gene masks, undefined constant profiles, explicit pair selection, retained flags, sorting alignment, and display-only jitter were checked. Five synthetic illustrative figures were added.
- Nine single-plate figures were added. The four QC panels were regenerated with identical QC values and neutral labels; five DE/enrichment panels retain their original pixels and numerical results. A rerun on the current development build changes some DE estimates, so the historical snapshots are not represented as exact current-version outputs.
- Fourteen bulk comparison figures were anonymized by editing vector text, preserving plotted values and geometry. The six historical evidence tables retain all 160 body cells after label substitution. Reference accession labels were removed from the public figures and prose.

The default bundled workflow and the detailed chapter examples passed. Enrichment with biological gene-set collections and STARsolo read-level QC require separate optional inputs. The small synthetic GMT checks validate the enrichment interface, not the historical pathway findings.

## Documentation checks

All twelve chapters passed browser checks at 1440 px and 390 px: no horizontal page overflow, chapter navigation present, and mobile navigation opens and closes. Desktop and mobile screenshots were visually reviewed. Every figure has a working link to its full-size asset. All 110 PNG files decode successfully; the static checker verifies their references.


The current static check passes 34 documents, 1,223 local links/assets, 116 Python-block syntax checks, and all 110 figure references. Run `python scripts/validate_site.py` after rebuilding. It checks local page/asset links, fragment targets, duplicate HTML IDs, runnable Python-block syntax, image descriptions, old sample identifiers, private filesystem references, and synthetic-data checksums. API signature displays are not executable examples.

The twelve Markdown chapters are the editable source for twelve generated HTML chapters. Six old HTML names redirect to the corresponding revised chapters. Figures and browser assets are local, and all internal site links are relative so the site works below a GitHub Pages repository path.

## Measured Drug-seq case

The added Drug-seq case uses a prespecified subset of three same-batch plates, three anonymized compounds at two doses, and all DMSO controls on those plates: 42 wells per plate, 126 total. It uses separate case-scoped Sample 1–3 labels; these are not the earlier worked or synthetic sample datasets.

- Every prepared value matches the source count after an explicit barcode join and shared-feature alignment. Physical coordinates were independently checked against the plate-layout API. Source identities, raw matrices, and private mapping configurations stay outside the public repository.
- Shared inputs contain 58,777 feature rows, including 86 ERCC spike features. Only the `ERCC-` spike prefix is removed; endogenous ERCC genes are preserved. A shared condition-aware filter retains 9,839 modeled features representing 9,835 unique symbols.
- All 126 wells and 72 DMSO controls are retained. One control per plate is flagged for review. QC metrics and the recurring low-depth control position were independently checked.
- The plate-plus-condition design has nine columns and rank nine. All 24 final contrast summaries (three separate plates plus the combined model, six conditions each) were recomputed from the gene-level results. Every model contains 9,839 tested and estimable features.
- Six independent first-plate refits reproduce the one-plate estimates. All 18 cross-plate effect correlations and significant-feature overlap results were independently verified. Counts are hashed and checked unchanged through QC, filtering of retained values, normalization, visualization correction, and DE.
- Fourteen new figures, eight summary CSVs, and sanitized analysis settings were inspected for numerical consistency, legibility, and source-identity removal. The chapter's Python summary-reading example executed successfully; all twelve pages passed browser checks at desktop and phone widths. The strongest combined contrast is explicitly an exploratory display choice.
- Preparation fixtures check shuffled metadata, nonshared features, count equality, incompatible conditions, fractional-count rejection, and mixed/missing cell models. Case analysis rejects incompatible units/exposure, incomplete designs, invalid coordinates, and unsupported exclusion IDs.

The three plates share a cell-preparation batch and repeat the same spatial layout. Treatment/dose and position are not fully separable; adding row to this selected fixed-effect design makes it rank deficient. The combined model includes the one-plate example, and the common gene filter was computed over all selected plates. These constraints are described in the chapter; no independent biological-batch or held-out-validation claim is made.

## Reproduction boundaries

The worked figures and tables are fixed analysis examples with Sample 1–4 labels. The bundled practice data were generated independently and do not reproduce those numerical values. The practice runner deliberately uses a smaller gene-filtering rule, explicit species review settings, and no automatic DMSO exclusions. Bulk practice outputs are descriptive and do not recreate the historical bootstraps or pseudogene audit. Chapter text explains these differences.

Some report functions warn when optional housekeeping genes are not present in the synthetic feature list; the warning does not invalidate the core QC output. Matching the displayed package version alone is not enough to reproduce an analysis. Record the exact compatible source revision, environment, inputs, parameters, exclusions, and gene-set versions for real work.


## Terminology revision — 27 October 2026

The tutorials distinguish within-plate repeatability, between-plate consistency within one batch, and cross-plate treatment-effect concordance. The measured three-plate case does not claim experimental reproducibility under changed conditions. The input guide explains experimental and computational reproducibility separately, including the limited role of correlations.

Public drugseqpy API and signature names are unchanged. Tutorial result variables and newly generated output paths use concordance terminology; legacy chapter anchors remain available. This editorial revision does not change analysis inputs, numerical tables, filters, fitted models, or existing plotted values.

Revision checks passed: the rebuilt site resolves all local links; all 116 Python examples parse; every CSV and PNG matches its pre-revision bytes. A normalized-AST comparison confirms that the two edited analysis runners retain the same calculations, changing only local names, the concordance output filename, and terminology metadata.

The synthetic comparison workflow completed successfully after the terminology revision (two successful contrasts) and produced `cross_plate_concordance.png`. All twelve chapters passed desktop/mobile navigation and overflow checks; existing chapter anchors remain valid.
