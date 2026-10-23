# Synthetic practice inputs

All counts and identifiers in these four bundles were generated independently by `scripts/make_practice_data.py`. They contain no experimental measurements. The seed and file checksums are in `practice-manifest.json`.

Each bundle contains 48 wells, 800 features, a sparse Matrix Market count matrix, feature and barcode lists, and metadata. Sample 1–2 demonstrate treatment comparison; Sample 3–4 demonstrate human/mouse species review. Synthetic feature IDs use recognizable prefixes only and are not a reference annotation. Marker symbols are illustrative.

Run the tutorials from the repository root. See [Getting started](../docs/getting-started.md). To create a new copy, use `python scripts/make_practice_data.py --output practice-copy`; the generator refuses to overwrite existing sample folders.
