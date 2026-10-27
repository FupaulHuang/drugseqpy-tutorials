#!/usr/bin/env python3
"""Run the connected tutorials with your own four anonymized input bundles.

Each data/sample_N folder contains metadata.csv and matrix/{matrix.mtx.gz,
features.tsv.gz,barcodes.tsv.gz}. Metadata must include barcode, well_id,
compound, dose and sample_type. Folder numbers define the displayed plate IDs.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import drugseqpy as ds


def load_samples(data_dir: Path, numbers: tuple[int, ...]):
    """Join on barcodes, assign public labels, and merge identical references."""
    if not numbers:
        raise ValueError("Select at least one input sample")
    objects = []
    for number in numbers:
        folder = data_dir / f"sample_{number}"
        metadata = pd.read_csv(folder / "metadata.csv")
        required = {"barcode", "well_id", "compound", "dose", "sample_type"}
        missing = sorted(required - set(metadata.columns))
        if missing:
            raise ValueError(f"{folder.name}/metadata.csv needs columns: {missing}")
        for column in ("barcode", "compound", "sample_type"):
            if metadata[column].isna().any() or metadata[column].astype(str).str.strip().eq("").any():
                raise ValueError(f"{folder.name}: {column} must be present for every well")
        if metadata["barcode"].duplicated().any():
            raise ValueError(f"{folder.name}: barcode must be unique within the plate")
        obj = ds.create_drugseq_from_files(
            folder / "matrix", folder / "metadata.csv",
            metadata_key="barcode", min_overlap=1.0,
        )
        label = f"Sample {number}"
        obj.obs["plate_id"] = label
        obj.obs["sample_id"] = label + "_" + obj.obs["barcode"].astype(str)
        obj.adata.obs_names = obj.obs["sample_id"].astype(str)
        if obj.obs["well_id"].isna().any() or obj.obs["well_id"].duplicated().any():
            raise ValueError(f"{label}: well_id must be present and unique within the plate")
        values = obj.layers["counts"].data
        if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
            raise ValueError(f"{label}: provide finite, nonnegative integer counts")
        if not obj.n_obs or not obj.n_vars or not (values > 0).any():
            raise ValueError(f"{label}: the input needs wells, features and at least one positive count")
        if objects and not obj.var_names.equals(objects[0].var_names):
            raise ValueError("Use the same ordered reference features for all plates in a workflow")
        if not obj.var_names.is_unique:
            raise ValueError(f"{label}: feature IDs must be unique")
        if objects and "gene_symbol" in obj.var and not obj.var["gene_symbol"].equals(objects[0].var["gene_symbol"]):
            raise ValueError("Gene symbols disagree between the input references")
        objects.append(obj)
    if len(objects) == 1:
        return objects[0]
    return ds.merge_drugseq_objects(*objects, batch_key="source_plate")


def save_figure(figure, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(figure)


def quality_control(dsd, output: Path, plate_rows: int, plate_columns: int):
    output.mkdir(parents=True, exist_ok=True)
    ds.compute_qc_metrics(dsd, inplace=True)
    audit = ds.validate_metadata(dsd, verbose=False)
    (output / "metadata_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    dsd.obs.to_csv(output / "sample_qc.csv", index=False)
    ds.generate_qc_metric_violin_figure(
        dsd, output_dir=output, formats=("png",),
        metrics=("total_umi", "n_genes_det", "pct_mito", "pct_ribo"),
        groupby="plate_id",
    )
    ds.generate_outlier_score_scatter_figure(
        dsd, output_dir=output, formats=("png",), groupby="plate_id",
        threshold=None, label_outliers=True,
    )
    ds.generate_plate_layout_figures(
        dsd, output_dir=output, formats=("png",), plate_value="total_umi",
        continuous=True, log_values=True, nrow=plate_rows, ncol=plate_columns,
        nrow_plate=1, output_name="plate_total_umi",
    )
    plt.close("all")


def species_workflow(dsd, output: Path, species_map: dict[str, str], threshold: float, minimum_umi: int):
    output.mkdir(parents=True, exist_ok=True)
    for prefix in ("ENSG", "ENSMUSG"):
        if not dsd.var_names.str.startswith(prefix).any():
            raise ValueError(f"No feature IDs start with {prefix}; use a joint human/mouse reference")
    ds.assign_well_species(dsd, source_col="compound", mapping=species_map)
    unknown = dsd.obs.loc[dsd.obs["expected_species"].eq("unknown"), "compound"].unique()
    if len(unknown):
        raise ValueError(f"Unmapped compound labels: {unknown.tolist()}; provide --species-map JSON")
    settings = dict(
        expected_col="expected_species",
        species_a_prefix="ENSG", species_b_prefix="ENSMUSG",
        species_a_name="human", species_b_name="mouse",
        assignment_threshold=threshold, min_species_umi=minimum_umi,
    )
    # The threshold and minimum depth are illustrative review settings.
    metrics = ds.compute_mix_species_qc(dsd, **settings)
    metrics.to_csv(output / "species_metrics.csv", index=False)
    ds.generate_mix_species_report(
        dsd, output_dir=output, formats=("png",), nrow_plate=2, **settings,
    )
    plt.close("all")


def apply_reviewed_exclusions(dsd, output: Path, exclusions: pd.DataFrame | None):
    """Remove only explicitly reviewed IDs, preserving the original QC audit."""
    if exclusions is not None:
        requested = exclusions["sample_id"].astype(str)
        unknown = sorted(set(requested) - set(dsd.obs["sample_id"].astype(str)))
        if unknown:
            raise ValueError(f"Exclusion IDs not found in this comparison: {unknown}")
        excluded = dsd.obs["sample_id"].astype(str).isin(requested)
        audit = dsd.obs.loc[excluded].reset_index(drop=True).merge(exclusions, on="sample_id", how="left", suffixes=("", "_review"))
        audit.to_csv(output / "reviewed_exclusions.csv", index=False)
        dsd = dsd.subset(obs_mask=(~excluded).to_numpy())
    return dsd


def check_comparison_result(result):
    manifest = result["manifest"]
    failures = manifest.loc[manifest["status"].isin(["failed", "partial"])]
    if not failures.empty or result["summary"].empty:
        raise RuntimeError("Comparison was incomplete; inspect the saved DE manifest")


def validate_comparison_conditions(analysis):
    """Each compound label must identify one complete dose/time condition."""
    for column in ("dose", "time"):
        if column not in analysis.obs:
            continue
        values = analysis.obs[column]
        if values.isna().any() or values.astype(str).str.strip().eq("").any():
            raise ValueError(f"Each biological well needs a complete {column} value")
        if analysis.obs.groupby("compound", observed=True)[column].nunique().gt(1).any():
            raise ValueError(
                f"Use one matching {column} value per compound across all selected wells and plates; "
                "give distinct conditions distinct compound labels"
            )


def single_plate_workflow(
    dsd, output: Path, full_embeddings: bool, gene_sets: dict[str, str],
    exclusions: pd.DataFrame | None, plate_rows: int, plate_columns: int,
):
    """Compare treatments with DMSO on one plate using the raw count layer."""
    if dsd.obs["plate_id"].nunique() != 1:
        raise ValueError("The single-plate workflow requires exactly one plate")
    quality_control(dsd, output / "qc", plate_rows, plate_columns)
    controls = ds.audit_control_outliers(dsd, control_label="DMSO", threshold=5.0)
    # Flags invite review; they are never automatic exclusion decisions.
    controls.rename(columns={"excluded": "flagged_for_review"}).to_csv(
        output / "control_review.csv", index=False,
    )
    n_input = dsd.n_obs
    dsd = apply_reviewed_exclusions(dsd, output, exclusions)
    excluded_labels = {"water", "blank", "media"}
    biological = ~dsd.obs["compound"].astype(str).str.casefold().isin(excluded_labels)
    analysis = dsd.subset(obs_mask=biological.to_numpy())
    design = pd.crosstab(analysis.obs["compound"], analysis.obs["plate_id"])
    design.to_csv(output / "design_counts.csv")
    if "DMSO" not in design.index or len(design.index) < 2:
        raise ValueError("Single-plate comparison needs compound='DMSO' and at least one treatment")
    if (design < 2).any().any():
        raise ValueError("Single-plate comparison needs >=2 wells per treatment and DMSO; inspect design_counts.csv")
    if (np.asarray(analysis.layers["counts"].sum(axis=1)).ravel() <= 0).any():
        raise ValueError("Comparison contains zero-count wells; review them before normalization and DE")
    validate_comparison_conditions(analysis)
    n_before_filter = analysis.n_vars
    # A deliberately modest rule for the small synthetic practice matrix.
    analysis = ds.filter_genes(analysis, min_count=1, min_samples=2, group_aware=True)
    if analysis.n_vars < 3:
        raise ValueError("Too few expressed genes remain for comparison")
    if (np.asarray(analysis.layers["counts"].sum(axis=1)).ravel() <= 0).any():
        raise ValueError("Gene filtering left a zero-count well; inspect the expression and design")
    analysis.obs.to_csv(output / "analysis_samples.csv", index=False)
    analysis.var.to_csv(output / "retained_features.csv", index_label="feature_id")
    (output / "processing_summary.json").write_text(json.dumps({
        "plate_id": str(analysis.obs["plate_id"].iloc[0]),
        "input_wells": n_input,
        "reviewed_exclusions": n_input - dsd.n_obs,
        "nonbiological_wells_removed": int((~biological).sum()),
        "analysis_wells": analysis.n_obs,
        "features_before_filter": n_before_filter,
        "features_after_filter": analysis.n_vars,
        "gene_filter": {"min_count": 1, "min_samples": 2, "group_aware": True},
        "normalization": "limma_voom",
        "batch_cols": None,
        "within_plate": False,
    }, indent=2) + "\n")
    # Normalization replaces X; the count layer remains the input to DE.
    ds.normalize_counts(analysis, method="limma_voom", inplace=True)
    n_pcs = min(20, analysis.n_obs - 1, analysis.n_vars - 1)
    ds.run_pca(analysis, n_pcs=n_pcs, use_highly_variable=False)
    save_figure(ds.plot_pca(analysis, color_by="compound"), output / "pca_compound.png")
    if full_embeddings:
        ds.run_umap(analysis, dims=n_pcs, n_neighbors=min(15, analysis.n_obs - 1), random_state=42)
        ds.run_tsne(analysis, dims=n_pcs, perplexity=min(30, (analysis.n_obs - 1) / 3), random_state=42)
        save_figure(ds.plot_umap(analysis, color_by="compound"), output / "umap_compound.png")
        save_figure(ds.plot_tsne(analysis, color_by="compound"), output / "tsne_compound.png")
    result = ds.run_comparison_workflow(
        analysis, output_dir=output / "de", group_col="compound", reference="DMSO",
        exclude_levels=("water", "blank", "media"),
        methods=("limma_voom",), batch_cols=None, within_plate=False,
        min_replicates=2, n_jobs=1,
        enrichment=bool(gene_sets), libraries=tuple(gene_sets), gene_set_paths=gene_sets,
        enrichment_mode=("ora",), background="all_tested", allow_network=False,
        figures=("volcano", "heatmap", "bar", "enrichment") if gene_sets else ("volcano", "heatmap", "bar"),
        n_top_genes=30, formats=("png",),
    )
    check_comparison_result(result)
    plt.close("all")


def comparison_workflow(dsd, output: Path, full_embeddings: bool, gene_sets: dict[str, str], exclusions: pd.DataFrame | None = None):
    output.mkdir(parents=True, exist_ok=True)
    ds.compute_qc_metrics(dsd, inplace=True)
    controls = ds.audit_control_outliers(dsd, control_label="DMSO", threshold=5.0)
    # The API names this flag "excluded", but the audit does not remove wells.
    controls.rename(columns={"excluded": "flagged_for_review"}).to_csv(
        output / "control_review.csv", index=False,
    )
    dsd = apply_reviewed_exclusions(dsd, output, exclusions)
    excluded_labels = {"water", "blank", "media"}
    biological = ~dsd.obs["compound"].astype(str).str.casefold().isin(excluded_labels)
    analysis = dsd.subset(obs_mask=biological.to_numpy())
    if (np.asarray(analysis.layers["counts"].sum(axis=1)).ravel() <= 0).any():
        raise ValueError("Comparison contains zero-count wells; review them before normalization and DE")
    design = pd.crosstab(analysis.obs["compound"], analysis.obs["plate_id"])
    design.to_csv(output / "design_counts.csv")
    if "DMSO" not in design.index or len(design.index) < 2:
        raise ValueError("Comparison needs compound='DMSO' and at least one treatment")
    if (design < 2).any().any():
        raise ValueError("This worked comparison needs >=2 wells per treatment and DMSO on each plate; inspect design_counts.csv")
    validate_comparison_conditions(analysis)
    analysis = ds.filter_genes(analysis, min_count=1, min_samples=2, group_aware=True)
    if analysis.n_vars < 3:
        raise ValueError("Too few expressed genes remain for comparison")

    corrected = analysis.copy()
    ds.correct_batch(corrected, batch_col="plate_id", method="limma_voom", covariate_cols=["compound"])
    before = corrected.copy()
    # Match the transformation exactly; compare only the effect of batch removal.
    before.adata.X = corrected.layers["batch_uncorrected"].copy()
    n_pcs = min(20, analysis.n_obs - 1, analysis.n_vars - 1)
    if n_pcs < 2:
        raise ValueError("At least three wells and genes are needed for these views")
    for name, obj in (("before", before), ("after", corrected)):
        ds.run_pca(obj, n_pcs=n_pcs, use_highly_variable=False)
        save_figure(ds.plot_pca(obj, color_by="plate_id"), output / f"pca_{name}_plate.png")
        save_figure(ds.plot_pca(obj, color_by="compound"), output / f"pca_{name}_compound.png")
        if full_embeddings:
            ds.run_umap(obj, dims=n_pcs, n_neighbors=min(15, obj.n_obs - 1), random_state=42)
            ds.run_tsne(obj, dims=n_pcs, perplexity=min(30, (obj.n_obs - 1) / 3), random_state=42)
            save_figure(ds.plot_umap(obj, color_by="plate_id"), output / f"umap_{name}.png")
            save_figure(ds.plot_tsne(obj, color_by="plate_id"), output / f"tsne_{name}.png")
    ds.compute_batch_effect_metrics(before, corrected).to_csv(output / "batch_metrics.csv", index=False)

    repeatability = ds.compute_within_plate_repeatability(analysis, use_norm=False)
    repeatability.to_csv(output / "within_plate_pairs.csv", index=False)
    ds.summarise_within_plate_repeatability(repeatability).to_csv(output / "within_plate_summary.csv", index=False)
    save_figure(ds.plot_within_plate_repeatability(repeatability), output / "within_plate_repeatability.png")

    # DE uses the raw count layer with plate included in the model.
    result = ds.run_comparison_workflow(
        analysis, output_dir=output / "de", group_col="compound", reference="DMSO",
        methods=("limma_voom",), batch_cols=["plate_id"], within_plate=True,
        min_replicates=2, n_jobs=1, normalize="limma_voom",
        enrichment=bool(gene_sets), libraries=tuple(gene_sets), gene_set_paths=gene_sets,
        enrichment_mode=("ora",), allow_network=False,
        figures=("volcano", "heatmap", "bar", "enrichment") if gene_sets else ("volcano", "heatmap", "bar"),
        n_top_genes=30, formats=("png",),
    )
    check_comparison_result(result)
    plt.close("all")
    # Public API name retained; results describe cross-plate effect concordance.
    concordance = ds.compute_cross_plate_reproducibility(analysis, exclude_levels=())
    concordance["summary"].to_csv(output / "cross_plate_summary.csv", index=False)
    for (plate, compound), table in concordance["plate_results"].items():
        filename = re.sub(r"[^A-Za-z0-9_-]+", "_", f"{plate}_{compound}")
        table.to_csv(output / f"{filename}_de.csv", index=False)
    save_figure(ds.plot_cross_plate_reproducibility(concordance["summary"]), output / "cross_plate_concordance.png")
    save_figure(ds.plot_plate_logfc_concordance(concordance), output / "cross_plate_logfc.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workflow", choices=("qc", "species", "single-plate", "comparison", "all"))
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--plate-rows", type=int, default=8)
    parser.add_argument("--plate-columns", type=int, default=12)
    parser.add_argument("--full-embeddings", action="store_true", help="Also run UMAP and t-SNE")
    parser.add_argument("--exclude-samples", type=Path, help="Reviewed exclusions for single-plate/comparison/all: CSV with sample_id and optional reason")
    parser.add_argument("--species-map", type=Path, help="JSON mapping compound labels to human, mouse or control")
    parser.add_argument("--species-threshold", type=float, default=0.9, help="Illustrative species-assignment fraction (default: 0.9)")
    parser.add_argument("--min-species-umi", type=int, default=100, help="Illustrative minimum assigned-species UMIs (default: 100)")
    parser.add_argument("--gene-set", action="append", default=[], metavar="LIBRARY=GMT", help="Optional local GMT, e.g. hallmark=gene_sets/hallmark.gmt")
    args = parser.parse_args()
    if not 0.5 < args.species_threshold <= 1.0:
        parser.error("--species-threshold must be greater than 0.5 and at most 1")
    if args.min_species_umi < 0:
        parser.error("--min-species-umi must be nonnegative")
    if args.plate_rows < 1 or args.plate_columns < 1:
        parser.error("Plate dimensions must be positive")
    exclusions = None
    if args.exclude_samples:
        if args.workflow not in {"single-plate", "comparison", "all"}:
            parser.error("--exclude-samples applies to single-plate, comparison or all")
        exclusions = pd.read_csv(args.exclude_samples, dtype={"sample_id": str})
        if "sample_id" not in exclusions or exclusions["sample_id"].isna().any() or exclusions["sample_id"].str.strip().eq("").any() or exclusions["sample_id"].duplicated().any():
            parser.error("Exclusion CSV requires unique, non-missing sample_id values")
    gene_sets = {}
    for entry in args.gene_set:
        name, separator, filename = entry.partition("=")
        if not separator or not name or not Path(filename).is_file():
            parser.error("--gene-set requires LIBRARY=existing-local-file.gmt")
        gene_sets[name] = filename
    species_map = {"231": "human", "3T3": "mouse", "water": "control", "blank": "control", "human": "human", "mouse": "mouse"}
    if args.species_map:
        species_map = json.loads(args.species_map.read_text())
    screen = None
    if args.workflow in {"single-plate", "comparison", "all"}:
        screen = load_samples(args.data_dir, (1,) if args.workflow == "single-plate" else (1, 2))
        if exclusions is not None:
            unknown = sorted(set(exclusions["sample_id"]) - set(screen.obs["sample_id"].astype(str)))
            if unknown:
                parser.error(f"Exclusion IDs not found in the selected comparison inputs: {unknown}")
    if args.workflow in {"qc", "species", "all"}:
        mixed = load_samples(args.data_dir, (3, 4))
        if args.workflow in {"qc", "all"}:
            quality_control(mixed, args.output_dir / "qc", args.plate_rows, args.plate_columns)
        if args.workflow in {"species", "all"}:
            species_workflow(mixed, args.output_dir / "species", species_map, args.species_threshold, args.min_species_umi)
    if args.workflow in {"single-plate", "all"}:
        single = screen.subset(obs_mask=screen.obs["plate_id"].eq("Sample 1").to_numpy())
        single_exclusions = exclusions
        if exclusions is not None:
            single_exclusions = exclusions.loc[exclusions["sample_id"].isin(single.obs["sample_id"])].copy()
        single_plate_workflow(
            single, args.output_dir / "single-plate", args.full_embeddings, gene_sets,
            single_exclusions, args.plate_rows, args.plate_columns,
        )
    if args.workflow in {"comparison", "all"}:
        comparison_workflow(screen, args.output_dir / "comparison", args.full_embeddings, gene_sets, exclusions)
    print(f"Completed {args.workflow}. Outputs: {args.output_dir}")


if __name__ == "__main__":
    main()
