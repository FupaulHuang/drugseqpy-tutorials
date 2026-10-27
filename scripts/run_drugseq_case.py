#!/usr/bin/env python3
"""Analyze anonymized, real UMI plates with matched treatment/dose conditions.

Inputs follow the tutorial sample_N/matrix + metadata.csv contract and must
also include condition, time and physical x/y coordinates. Outputs contain no raw count matrix.
No wells are removed unless explicitly listed in --exclude-samples.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from PIL import Image

import drugseqpy as ds
from run_tutorial import load_samples


def count_digest(obj):
    matrix = sparse.csr_matrix(obj.layers["counts"]).copy()
    matrix.sort_indices()
    digest = hashlib.sha256()
    for value in (np.asarray(matrix.shape, dtype="int64"), matrix.indptr,
                  matrix.indices, matrix.data):
        digest.update(value.tobytes())
    return digest.hexdigest()


def save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    # Strip PNG metadata from every shareable figure.
    with Image.open(path) as source:
        clean = Image.frombytes(source.mode, source.size, source.tobytes())
    clean.save(path, optimize=True)


def clean_generated_png(path):
    with Image.open(path) as source:
        clean = Image.frombytes(source.mode, source.size, source.tobytes())
    clean.save(path, optimize=True)


def check_design(obj):
    required = ["condition", "compound", "dose", "dose_unit", "time", "row", "col"]
    missing = set(required) - set(obj.obs)
    if missing or obj.obs[required].isna().any().any():
        raise ValueError(f"Complete condition/dose/time/physical-coordinate metadata required; missing: {missing}")
    for column in required:
        if obj.obs[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"Empty {column} values")
    if obj.obs.groupby("condition", observed=True)[["compound", "dose", "dose_unit", "time"]].nunique().gt(1).any().any():
        raise ValueError("Every condition must identify exactly one compound, dose and time")
    dose = pd.to_numeric(obj.obs.dose, errors="raise")
    if not np.isfinite(dose).all() or (dose < 0).any():
        raise ValueError("Doses must be finite nonnegative numbers")
    if obj.obs.time.nunique() != 1 or obj.obs.dose_unit.nunique() != 1:
        raise ValueError("This runner requires one shared exposure time and one dose unit, including DMSO")
    if not obj.obs.condition.eq("DMSO").equals(obj.obs.compound.eq("DMSO")):
        raise ValueError("DMSO controls must have both condition and compound set to DMSO")
    design = pd.crosstab(obj.obs.condition, obj.obs.plate_id)
    if "DMSO" not in design.index or len(design) < 2 or (design < 2).any().any():
        raise ValueError("Every selected plate needs DMSO and >=2 wells for every selected treatment condition")
    dummy = pd.get_dummies(obj.obs[["plate_id", "condition"]], drop_first=True, dtype=float)
    matrix = np.column_stack([np.ones(obj.n_obs), dummy.to_numpy()])
    rank = int(np.linalg.matrix_rank(matrix))
    if rank != matrix.shape[1]:
        raise ValueError("Plate plus condition model is rank deficient")
    return design, {"columns": ["intercept", *dummy.columns], "rank": rank,
                    "n_columns": int(matrix.shape[1]), "n_observations": obj.n_obs}


def pca_panels(before, after, color_by, path):
    groups = sorted(before.obs[color_by].astype(str).unique())
    colors = plt.get_cmap("tab10")(np.linspace(0, .9, len(groups)))
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    points = [np.asarray(x.adata.obsm["X_pca"])[:, :2] for x in (before, after)]
    low = np.vstack(points).min(axis=0); high = np.vstack(points).max(axis=0)
    pad = np.maximum((high - low) * .08, .1)
    for ax, obj, coords, title in zip(axes, (before, after), points, ("Before plate adjustment", "After plate adjustment")):
        for group, color in zip(groups, colors):
            mask = obj.obs[color_by].astype(str).eq(group).to_numpy()
            ax.scatter(coords[mask, 0], coords[mask, 1], label=group, c=[color], s=32, alpha=.85)
        ax.set(title=title, xlabel="PC1", ylabel="PC2",
               xlim=(low[0]-pad[0], high[0]+pad[0]), ylim=(low[1]-pad[1], high[1]+pad[1]))
        ax.grid(alpha=.18)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=min(4, len(groups)), frameon=False)
    fig.suptitle("Same log-expression transformation; PCA fitted separately for each view")
    fig.tight_layout(rect=(0, .12 if len(groups) > 4 else .07, 1, .94))
    save(fig, path)


def de(obj, *, batch_cols, output, label):
    ds.compute_multi_de(obj, group_col="condition", reference="DMSO",
                        methods=("limma_voom",), batch_cols=batch_cols,
                        within_plate=bool(batch_cols), exclude_levels=(),
                        min_replicates=2, n_jobs=1, fdr_threshold=.05,
                        lfc_threshold=.5, inplace=True)
    failures = obj.adata.uns.get("comparison_failures", {})
    tables = obj.adata.uns.get("comparison_results", {})
    expected = set(obj.obs.condition) - {"DMSO"}
    if failures or set(tables) != expected:
        raise RuntimeError(f"Incomplete {label} DE: {failures}")
    output.mkdir(parents=True, exist_ok=True)
    for condition, table in tables.items():
        table.to_csv(output / f"{condition}.csv", index=False)
    summary = ds.summarise_comparison(obj)
    summary["n_estimable"] = summary.contrast.map({
        name: int((table.padj.notna() & table.logFC.notna()).sum()) for name, table in tables.items()
    })
    summary.insert(0, "analysis", label)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--sample-numbers", type=int, nargs="+", default=[1, 2, 3])
    parser.add_argument("--output-dir", type=Path, default=Path("results/drugseq-case"))
    parser.add_argument("--exclude-samples", type=Path, help="Reviewed CSV: unique sample_id, optional reason")
    args = parser.parse_args()
    if len(args.sample_numbers) < 2 or len(set(args.sample_numbers)) != len(args.sample_numbers) or min(args.sample_numbers) < 1:
        parser.error("Choose at least two distinct sample numbers")
    output = args.output_dir
    figures = output / "figures"; figures.mkdir(parents=True, exist_ok=True)
    tables_dir = output / "tables"; tables_dir.mkdir(parents=True, exist_ok=True)
    obj = load_samples(args.data_dir, tuple(args.sample_numbers))
    if not {"x", "y"}.issubset(obj.obs):
        raise ValueError("Explicit physical x (column) and y (row) coordinates are required")
    obj.obs["row"] = pd.to_numeric(obj.obs["y"], errors="raise")
    obj.obs["col"] = pd.to_numeric(obj.obs["x"], errors="raise")
    if not (obj.obs.row.between(1, 16) & obj.obs.col.between(1, 24)).all():
        raise ValueError("Coordinates must fit the physical 16 by 24 plate")
    if not (obj.obs.row.mod(1).eq(0) & obj.obs.col.mod(1).eq(0)).all():
        raise ValueError("Physical plate coordinates must be integers")
    if obj.obs.duplicated(["plate_id", "row", "col"]).any():
        raise ValueError("Physical well positions must be unique within a plate")
    check_design(obj)
    initial_digest = count_digest(obj)
    input_wells, input_features = obj.n_obs, obj.n_vars
    symbols = obj.var["gene_symbol"].astype(str) if "gene_symbol" in obj.var else pd.Series(obj.var_names, index=obj.var_names)
    spike = symbols.str.startswith("ERCC-").to_numpy()
    ds.compute_qc_metrics(obj, inplace=True)
    obj.obs["biological_umi"] = np.asarray(obj.layers["counts"][:, ~spike].sum(axis=1)).ravel()
    obj.obs["ercc_umi"] = np.asarray(obj.layers["counts"][:, spike].sum(axis=1)).ravel()
    obj.obs["biological_genes_detected"] = np.asarray((obj.layers["counts"][:, ~spike] > 0).sum(axis=1)).ravel()
    assert count_digest(obj) == initial_digest
    obj.obs.to_csv(output / "selected_well_qc.csv", index=False)
    qc = obj.obs.groupby("plate_id", observed=True).agg(
        wells=("sample_id", "size"), median_total_umi=("total_umi", "median"),
        minimum_total_umi=("total_umi", "min"), median_biological_umi=("biological_umi", "median"),
        median_detected_features=("n_genes_det", "median"), median_mito_pct=("pct_mito", "median"),
        median_rpl_rps_pct=("pct_ribo", "median"),
    ).reset_index()
    qc.to_csv(tables_dir / "qc-summary.csv", index=False)
    ds.generate_qc_metric_violin_figure(obj, output_dir=figures, formats=("png",),
        report_label="Selected_Plates", metrics=("total_umi", "n_genes_det", "pct_mito", "pct_ribo"), groupby="plate_id")
    # The package names multi-plate output consistently; discover by its suffix.
    candidates = list(figures.glob("*QC_Boxplots.png"))
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one QC figure, found {candidates}")
    candidates[0].replace(figures / "drugseq-case-01.png")
    for number, value, continuous in [(2, "condition", False), (3, "biological_umi", True)]:
        ds.generate_plate_layout_figures(obj, output_dir=figures, formats=("png",),
            plate_value=value, continuous=continuous, nrow=16, ncol=24, nrow_plate=3,
            output_name=f"drugseq-case-{number:02}", figsize=(25, 7.8), log_values=True)
    control = ds.audit_control_outliers(obj, control_col="condition", threshold=5.0)
    control = control.rename(columns={"excluded": "flagged_for_review", "exclusion_rule": "review_rule"})
    control.to_csv(output / "control_review.csv", index=False)
    review_counts = control.groupby("plate_id", observed=True).agg(
        controls=("sample_id", "size"), flagged_for_review=("flagged_for_review", "sum"))
    review_counts.to_csv(tables_dir / "control-review-summary.csv")
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for plate, group in control.groupby("plate_id", observed=True):
        ax.scatter(group.total_umi, group.n_genes_det, label=plate, s=35, alpha=.8)
    flagged = control[control.flagged_for_review]
    ax.scatter(flagged.total_umi, flagged.n_genes_det, marker="o", facecolors="none", edgecolors="#b4233b", s=110, label="Flagged for review; retained")
    ax.set(xlabel="Total UMI, including ERCC spikes", ylabel="Detected features, including ERCC spikes",
           title="DMSO control review: no automatic exclusions")
    ax.legend(frameon=False); ax.grid(alpha=.18); fig.tight_layout(); save(fig, figures / "drugseq-case-04.png")
    excluded_count = 0
    if args.exclude_samples:
        exclusions = pd.read_csv(args.exclude_samples, dtype={"sample_id": str})
        if ("sample_id" not in exclusions or exclusions.sample_id.isna().any()
                or exclusions.sample_id.str.strip().eq("").any() or exclusions.sample_id.duplicated().any()):
            raise ValueError("Exclusions require unique, nonmissing sample_id values")
        unknown = set(exclusions.sample_id) - set(obj.obs.sample_id)
        if unknown:
            raise ValueError(f"Exclusions not found: {sorted(unknown)}")
        selected = obj.obs.sample_id.isin(exclusions.sample_id)
        excluded_count = int(selected.sum())
        obj.obs.loc[selected].reset_index(drop=True).merge(exclusions, on="sample_id", suffixes=("", "_review")).to_csv(output / "reviewed_exclusions.csv", index=False)
        obj = obj.subset(obs_mask=(~selected).to_numpy())
    design, rank = check_design(obj)
    design.to_csv(tables_dir / "design-counts.csv")
    spatial = obj.obs.groupby(["plate_id", "condition", "row"], observed=True).agg(wells=("sample_id", "size")).reset_index()
    spatial.to_csv(tables_dir / "spatial-design.csv", index=False)
    biological = obj.subset(var_mask=~spike)
    analysis = ds.filter_genes(biological, min_count=5, min_samples=3,
        group_aware=True, group_col="condition")
    expected_counts = biological.layers["counts"][:, biological.var_names.get_indexer(analysis.var_names)]
    assert (sparse.csr_matrix(expected_counts) != sparse.csr_matrix(analysis.layers["counts"])).nnz == 0
    if (np.asarray(analysis.layers["counts"].sum(axis=1)).ravel() <= 0).any():
        raise ValueError("A selected well has no retained biological counts; review it explicitly")
    analysis_digest = count_digest(analysis)
    analysis.var.to_csv(output / "retained_features.csv", index_label="feature_id")
    before = analysis.copy(); corrected = analysis.copy()
    ds.correct_batch(corrected, batch_col="plate_id", covariate_cols=["condition"], method="limma_voom", inplace=True)
    before.adata.X = corrected.layers["batch_uncorrected"].copy()
    n_pcs = min(20, analysis.n_obs-1, analysis.n_vars-1)
    for view in (before, corrected):
        ds.run_pca(view, n_pcs=n_pcs, use_highly_variable=False, inplace=True)
        assert count_digest(view) == analysis_digest
    pca_panels(before, corrected, "plate_id", figures / "drugseq-case-05.png")
    pca_panels(before, corrected, "condition", figures / "drugseq-case-06.png")
    ds.compute_batch_effect_metrics(before, corrected, batch_col="plate_id", group_col="condition").to_csv(tables_dir / "batch-diagnostics.csv", index=False)
    repeat = ds.compute_within_plate_repeatability(analysis, group_col="condition", use_norm=False)
    repeat.to_csv(output / "repeatability_pairs.csv", index=False)
    ds.summarise_within_plate_repeatability(repeat, group_col="condition").to_csv(tables_dir / "repeatability-summary.csv", index=False)
    repeat_fig = ds.plot_within_plate_repeatability(repeat, group_col="condition")
    repeat_fig.axes[0].set(title="Within-plate profile correlation by treatment-dose condition",
                          xlabel="Treatment-dose condition")
    save(repeat_fig, figures / "drugseq-case-07.png")
    first_plate = f"Sample {args.sample_numbers[0]}"
    chosen_condition = sorted(set(analysis.obs.condition) - {"DMSO"})[0]
    save(ds.plot_within_plate_correlation_heatmap(analysis, plate=first_plate, compound=chosen_condition,
        group_col="condition", use_norm=False), figures / "drugseq-case-08.png")
    save(ds.plot_within_plate_replicate_scatter(analysis, plate=first_plate, compound=chosen_condition,
        group_col="condition", use_norm=False, max_pairs=3), figures / "drugseq-case-09.png")
    single = analysis.subset(obs_mask=analysis.obs.plate_id.eq(first_plate).to_numpy())
    single_digest = count_digest(single)
    single_summary = de(single, batch_cols=None, output=output / "de" / "single", label=first_plate)
    combined_summary = de(analysis, batch_cols=["plate_id"], output=output / "de" / "combined", label="Combined plates")
    assert count_digest(single) == single_digest and count_digest(analysis) == analysis_digest
    # Independent per-plate DE uses the same retained feature universe and raw counts.
    ds.normalize_counts(analysis, method="limma_voom", inplace=True)
    # Public API name retained; results describe cross-plate effect concordance.
    concordance = ds.compute_cross_plate_reproducibility(analysis, group_col="condition", reference="DMSO",
        exclude_levels=(), method="limma_voom", min_replicates=2, fdr_threshold=.05, lfc_threshold=.5)
    assert count_digest(analysis) == analysis_digest
    concordance["summary"].to_csv(tables_dir / "cross-plate-summary.csv", index=False)
    per_plate = []
    for (plate, condition), table in concordance["plate_results"].items():
        path = output / "de" / re.sub(r"\W+", "_", plate); path.mkdir(parents=True, exist_ok=True)
        table.to_csv(path / f"{condition}.csv", index=False)
        significant = (table.padj < .05) & (table.logFC.abs() >= .5)
        per_plate.append({"analysis": plate, "contrast": condition, "n_tested": len(table),
            "n_estimable": int((table.padj.notna() & table.logFC.notna()).sum()),
            "n_sig_up": int((significant & table.logFC.gt(0)).sum()),
            "n_sig_down": int((significant & table.logFC.lt(0)).sum()),
            "n_sig_total": int(significant.sum()), "min_padj": float(table.padj.min())})
        if plate == first_plate:
            original = single.adata.uns["comparison_results"][condition].set_index("gene")
            aligned = table.set_index("gene").loc[original.index]
            assert np.allclose(original.logFC, aligned.logFC, equal_nan=True)
            assert np.allclose(original.padj, aligned.padj, equal_nan=True)
    all_summary = pd.concat([pd.DataFrame(per_plate), combined_summary], ignore_index=True)
    all_summary.to_csv(tables_dir / "de-summary.csv", index=False)
    single_summary.to_csv(output / "single_plate_de_summary.csv", index=False)
    pivot = all_summary.pivot(index="contrast", columns="analysis", values="n_sig_total")
    fig, ax = plt.subplots(figsize=(12, 6))
    pivot.plot.bar(ax=ax); ax.set(ylabel="Features: FDR < 0.05 and |log2FC| >= 0.5", xlabel="Treatment-dose condition",
        title="Independent plate fits and combined plate-adjusted fit")
    ax.legend(frameon=False); ax.tick_params(axis="x", rotation=30); fig.tight_layout(); save(fig, figures / "drugseq-case-10.png")
    priority = str(combined_summary.sort_values(["n_sig_total", "min_padj", "contrast"], ascending=[False, True, True]).iloc[0].contrast)
    for number, fitted, title in [(11, analysis, "Combined plates"), (12, single, first_plate)]:
        fig = ds.plot_comparison_volcano(fitted, contrast=priority, n_label=0, figsize=(7.5, 6))
        fig.axes[0].set_title(f"{title}: {priority} vs DMSO\nExploratory contrast selected by combined DE feature count")
        save(fig, figures / f"drugseq-case-{number:02}.png")
    save(ds.plot_plate_logfc_concordance(concordance, compounds=[priority], max_panels=3,
        value_type="log2fc", figsize=(12, 10)), figures / "drugseq-case-13.png")
    cross = concordance["summary"].copy()
    cross["pair"] = cross.plate_a + " / " + cross.plate_b
    heat = cross.pivot(index="compound", columns="pair", values="pearson_logfc")
    fig, ax = plt.subplots(figsize=(9.5, 6.5)); im = ax.imshow(heat.to_numpy(), vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set(xticks=np.arange(len(heat.columns)), xticklabels=heat.columns,
        yticks=np.arange(len(heat)), yticklabels=heat.index,
        title="Independent plate-specific log2FC agreement across all conditions")
    ax.set_xticklabels(heat.columns, rotation=25, ha="right")
    for row in range(len(heat)):
        for col in range(len(heat.columns)):
            value = heat.iloc[row, col]
            ax.text(col, row, f"{value:.3f}", ha="center", va="center", color="white" if abs(value)>.65 else "black")
    fig.colorbar(im, ax=ax, label="Pearson r across retained features"); fig.tight_layout(); save(fig, figures / "drugseq-case-14.png")
    for path in figures.glob("drugseq-case-*.png"):
        clean_generated_png(path)
    settings = {
        "samples": [f"Sample {n}" for n in args.sample_numbers], "input_wells": input_wells,
        "retained_wells": analysis.n_obs, "reviewed_exclusions": excluded_count,
        "features_input": input_features, "spike_features_excluded": int(spike.sum()),
        "spike_rule": "gene_symbol starts ERCC-; endogenous ERCC genes retained",
        "biological_features_before_filter": biological.n_vars, "features_after_filter": analysis.n_vars,
        "gene_filter": {"min_count": 5, "min_samples": 3, "group_aware": True, "group_col": "condition"},
        "shared_filter_scope": "all selected plates; identical retained input genes in every DE fit",
        "de": {"method": "limma_voom", "reference": "DMSO", "fdr": .05, "absolute_log2fc": .5,
               "single_plate_covariates": [], "combined_covariates": ["plate_id"]},
        "design": rank, "exploratory_contrast": priority,
        "repeatability_illustration": {"plate": first_plate, "condition": chosen_condition,
            "selection": "first condition in lexical order; no correlation-based selection", "gene_universe": "all retained filtered genes, including shared zeros"},
        "counts_preserved": True, "filtered_counts_sha256": analysis_digest,
        "independent_first_plate_refit_agrees": True,
        "terminology": {
            "within_plate": "within-plate repeatability",
            "between_plate": "between-plate consistency within one cell-preparation batch",
            "effect_comparison": "cross-plate treatment-effect concordance",
            "experimental_reproducibility": "not established under changed experimental conditions",
        },
        "qc_scope": "selected wells only; total UMI/detected features include ERCC; biological_umi excludes ERCC-",
        "qc_annotation": "mitochondrial MT- and ribosomal-associated RPL/RPS symbol prefixes",
        "control_review_score": "outlier scores computed across all selected wells, including treatments; threshold 5 is only a review heuristic",
        "correction_scope": "visual diagnostics only; before view is exact batch_uncorrected transformation",
        "enrichment": "not run; no compatible local gene set supplied",
        "limitations": ["Selected plates share one cell preparation batch; no independent cell-batch validation",
            "Independent culture replication is unverified; wells are observations, not claimed biological replicates",
            "Combined analysis includes the first plate and is not external validation",
            "Repeated fixed layouts partially confound treatment/dose with position; plate adjustment does not remove this",
            "Strongest combined contrast is exploratory; six contrasts have separate within-contrast FDR adjustment"],
        "versions": {name: importlib.metadata.version(name) for name in ["drugseqpy", "numpy", "pandas", "scipy"]},
    }
    (output / "analysis_settings.json").write_text(json.dumps(settings, indent=2) + "\n")
    print(json.dumps({key: settings[key] for key in ["retained_wells", "features_after_filter", "exploratory_contrast", "counts_preserved"]}, indent=2))


if __name__ == "__main__":
    main()
