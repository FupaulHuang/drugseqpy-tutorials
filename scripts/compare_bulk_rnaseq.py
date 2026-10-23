#!/usr/bin/env python3
"""Descriptive, species-matched comparison of Sample 3/4 with bulk RNA-seq.

Bulk expression is a TSV with feature_id followed by reference columns.
Bulk metadata is a TSV with reference_id, study_id, species, compound, unit,
and primary columns. Blank expression cells mean unmeasured, never zero.
The bundled reference is entirely synthetic and has no biological meaning.
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
from scipy import sparse
from scipy.stats import spearmanr

import drugseqpy as ds
from run_tutorial import load_samples, save_figure


def stable_gene_ids(values, label):
    ids = pd.Index(values).astype(str).str.replace(r"\.\d+$", "", regex=True)
    if not ids.str.fullmatch(r"ENS(?:MUS)?G\d+").all():
        raise ValueError(f"{label}: use human/mouse Ensembl gene IDs, optionally with version suffixes")
    if not ids.is_unique:
        raise ValueError(f"{label}: duplicate stable gene IDs after stripping versions; resolve duplicates explicitly")
    return ids


def read_bulk(expression_path: Path, metadata_path: Path):
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str)
    required = {"reference_id", "study_id", "species", "compound", "unit", "primary"}
    missing = required - set(metadata)
    if missing:
        raise ValueError(f"Bulk metadata is missing columns: {sorted(missing)}")
    if metadata.empty or metadata[list(required)].isna().any().any():
        raise ValueError("Bulk metadata requires nonmissing values in every required column")
    if metadata[list(required)].apply(lambda s: s.str.strip().eq("").any()).any():
        raise ValueError("Bulk metadata required columns must not contain empty values")
    if metadata.reference_id.duplicated().any():
        raise ValueError("Bulk metadata reference_id must be unique")
    if not metadata.species.isin(["human", "mouse"]).all():
        raise ValueError("Bulk species must be human or mouse")
    if not metadata.unit.isin(["counts", "FPKM"]).all():
        raise ValueError("Bulk unit must be counts or FPKM")
    primary = metadata.primary.str.casefold().map({"true": True, "false": False, "1": True, "0": False})
    if primary.isna().any():
        raise ValueError("Bulk primary must be true/false or 1/0")
    metadata["primary"] = primary.astype(bool)
    if metadata.groupby("study_id").unit.nunique().gt(1).any():
        raise ValueError("Each study must use one expression unit")
    if metadata.groupby("compound").species.nunique().gt(1).any():
        raise ValueError("A compound label cannot refer to both species")
    expression = pd.read_csv(expression_path, sep="\t")
    if "feature_id" not in expression:
        raise ValueError("Bulk expression needs a feature_id column")
    expression.index = stable_gene_ids(expression.pop("feature_id"), "Bulk expression")
    if set(expression.columns) != set(metadata.reference_id):
        raise ValueError("Bulk expression columns must exactly match metadata reference_id values")
    expression = expression.apply(pd.to_numeric, errors="raise")
    values = expression.to_numpy(dtype=float)
    # Missing observations are allowed so a common measured universe can be formed.
    if np.isinf(values).any() or (values[np.isfinite(values)] < 0).any():
        raise ValueError("Observed bulk expression must be finite and nonnegative")
    for row in metadata.itertuples(index=False):
        observed = expression[row.reference_id].dropna().to_numpy()
        if not len(observed):
            raise ValueError(f"Reference {row.reference_id} has no measured expression")
        if row.unit == "counts" and not np.equal(observed, np.floor(observed)).all():
            raise ValueError(f"Reference {row.reference_id}: counts must be integers")
    return expression, metadata


def correlations(local, bulk):
    """Keep zeros and the fixed gene universe; a constant vector is undefined."""
    if np.ptp(local) == 0 or np.ptp(bulk) == 0:
        return {"spearman_rho": np.nan, "pearson_log1p": np.nan, "status": "constant_expression"}
    return {
        "spearman_rho": float(spearmanr(local, bulk).statistic),
        "pearson_log1p": float(np.corrcoef(np.log1p(local), np.log1p(bulk))[0, 1]),
        "status": "ok",
    }


def compare(data_dir: Path, expression_path: Path, metadata_path: Path, output: Path,
            species_map: dict[str, str], exclusions: pd.DataFrame | None = None):
    bulk, metadata = read_bulk(expression_path, metadata_path)
    plates = {f"Sample {number}": load_samples(data_dir, (number,)) for number in (3, 4)}
    for label, obj in plates.items():
        obj.adata.var_names = stable_gene_ids(obj.var_names, label)
        ds.assign_well_species(obj, source_col="compound", mapping=species_map)
    output.mkdir(parents=True, exist_ok=True)
    if exclusions is not None:
        known = set().union(*(set(obj.obs.sample_id) for obj in plates.values()))
        unknown = sorted(set(exclusions.sample_id) - known)
        if unknown:
            raise ValueError(f"Reviewed exclusions not found in Sample 3/4: {unknown}")
        exclusion_audits = []
        for label, obj in plates.items():
            selected = exclusions.loc[exclusions.sample_id.isin(obj.obs.sample_id)]
            mask = obj.obs.sample_id.isin(selected.sample_id)
            exclusion_audits.append(obj.obs.loc[mask].reset_index(drop=True).merge(selected, on="sample_id", suffixes=("", "_review")))
            plates[label] = obj.subset(obs_mask=(~mask).to_numpy())
        pd.concat(exclusion_audits, ignore_index=True).to_csv(output / "reviewed_exclusions.csv", index=False)
    profiles, wells, included, gene_sets = [], [], [], []
    for (species, compound), references in metadata.groupby(["species", "compound"], sort=True):
        if species_map.get(compound) != species:
            raise ValueError(f"Compound {compound!r} does not map to {species}; provide a matching --species-map JSON")
        matching = {
            label: obj.obs.loc[obj.obs.compound.astype(str).eq(compound) & obj.obs.expected_species.eq(species)]
            for label, obj in plates.items()
        }
        for label, observations in matching.items():
            if observations.empty:
                raise ValueError(f"{label}: no expected {species} wells with compound={compound!r}")
        # A compound label must identify one matched local condition, not a
        # pool of doses or time points that differ between the two plates.
        pooled_metadata = pd.concat(matching.values())
        for condition in ("dose", "time"):
            if condition in pooled_metadata:
                values = pooled_metadata[condition]
                if values.isna().any() or values.astype(str).str.strip().eq("").any() or values.nunique() != 1:
                    raise ValueError(
                        f"{species}/{compound}: use one complete, matching {condition} value across both samples; "
                        "encode distinct conditions as distinct compound labels"
                    )
        prefix = "ENSG" if species == "human" else "ENSMUSG"
        reference_ids = references.reference_id.tolist()
        shared = bulk.index[bulk.index.str.startswith(prefix) & bulk[reference_ids].notna().all(axis=1)]
        for obj in plates.values():
            shared = shared.intersection(obj.var_names, sort=False)
        shared = shared.sort_values()
        if len(shared) < 3:
            raise ValueError(f"{species}/{compound}: fewer than three genes measured in every plate/reference")
        gene_sets.extend({"species": species, "compound": compound, "feature_id": gene} for gene in shared)
        profile_values = {}
        for label, obj in plates.items():
            mask = obj.obs.compound.astype(str).eq(compound) & obj.obs.expected_species.eq(species)
            observations = matching[label]
            included.append(observations[["sample_id", "plate_id", "well_id", "compound", "expected_species"]])
            matrix = obj.layers["counts"][mask.to_numpy()][:, obj.var_names.get_indexer(shared)]
            matrix = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
            mean_counts = matrix.mean(axis=0)
            profile_values[label] = mean_counts
            for reference in references.itertuples(index=False):
                values = bulk.loc[shared, reference.reference_id].to_numpy(dtype=float)
                detail = {
                    "species": species, "compound": compound,
                    "reference_id": reference.reference_id, "study_id": reference.study_id,
                    "unit": reference.unit, "primary": reference.primary,
                    "plate_id": label, "n_genes": len(shared),
                }
                metric = correlations(mean_counts, values)
                if metric["status"] != "ok":
                    raise ValueError(f"{label}/{reference.reference_id}: mean profile or bulk reference is constant")
                profiles.append({**detail, "n_wells": len(observations), **metric})
                for row, (sample_id, well_id) in enumerate(zip(observations.sample_id, observations.well_id)):
                    wells.append({**detail, "sample_id": sample_id, "well_id": well_id,
                                  **correlations(matrix[row], values)})
        # Separate axis labels retain the original units; no common-unit claim.
        first = references.iloc[0]
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), sharex=True, sharey=True)
        x = bulk.loc[shared, first.reference_id].to_numpy(dtype=float)
        for ax, (label, values) in zip(axes, profile_values.items()):
            rho = correlations(values, x)["spearman_rho"]
            ax.scatter(np.log1p(x), np.log1p(values), s=11, alpha=.55, color="#147d75")
            ax.set(title=f"{label}; Spearman rho = {rho:.3f}",
                   xlabel=f"log1p(bulk {first.unit})", ylabel="log1p(mean UMI per retained well)")
            ax.grid(alpha=.15)
        fig.suptitle(f"{species} / {compound} / {first.reference_id}; {len(shared)} shared genes")
        fig.tight_layout()
        slug = re.sub(r"[^A-Za-z0-9_-]+", "_", f"{species}_{compound}")
        save_figure(fig, output / f"scatter_{slug}.png")
    profile_table = pd.DataFrame(profiles)
    profile_table.to_csv(output / "profile_correlations.csv", index=False)
    pd.DataFrame(wells).to_csv(output / "well_correlations.csv", index=False)
    pd.concat(included).drop_duplicates().to_csv(output / "included_wells.csv", index=False)
    pd.DataFrame(gene_sets).to_csv(output / "shared_features.csv", index=False)
    keys = ["species", "compound", "reference_id", "study_id", "unit", "primary", "n_genes"]
    paired = profile_table.pivot(index=keys, columns="plate_id", values="spearman_rho").reset_index()
    # Clip only for a finite descriptive Fisher transform, not for inference.
    paired["delta_fisher_z"] = (np.arctanh(paired["Sample 4"].clip(-.999999, .999999))
                                 - np.arctanh(paired["Sample 3"].clip(-.999999, .999999)))
    paired.to_csv(output / "paired_reference_differences.csv", index=False)
    study_keys = ["species", "compound", "primary", "study_id"]
    studies = paired.groupby(study_keys, as_index=False).agg(
        delta_fisher_z=("delta_fisher_z", "mean"), n_references=("reference_id", "size"),
    )
    studies.to_csv(output / "study_differences.csv", index=False)
    studies.groupby(["species", "compound", "primary"], as_index=False).agg(
        study_equal_delta_fisher_z=("delta_fisher_z", "mean"),
        n_studies=("study_id", "size"), n_references=("n_references", "sum"),
    ).to_csv(output / "study_equal_summary.csv", index=False)
    (output / "comparison_settings.json").write_text(json.dumps({
        "plates": ["Sample 3", "Sample 4"],
        "primary_metric": "Spearman correlation on a fixed shared gene universe",
        "secondary_metric": "Pearson correlation of separately log1p-transformed expression",
        "local_profile": "mean raw UMI counts over retained wells of matching expected species and compound",
        "feature_universe": "intersection measured in both plates and all references for each species/compound; zeros retained",
        "duplicate_stable_ids": "error; no automatic aggregation",
        "difference": "Fisher z(Sample 4 rho) minus Fisher z(Sample 3 rho)",
        "fisher_clip": .999999,
        "study_weighting": "equal mean of within-study reference differences; primary and sensitivity references separate",
        "interpretation": "descriptive only; wells and genes are not independent plate or method replicates",
        "condition_matching": "local compound, dose and supplied time labels must match; biological matching to bulk references is the reader's responsibility",
        "batch_correction": False,
    }, indent=2) + "\n")
    print(f"Compared {len(profile_table)} matched profiles. Outputs: {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--bulk-expression", type=Path, default=Path("data/bulk_practice/expression.tsv"))
    parser.add_argument("--bulk-metadata", type=Path, default=Path("data/bulk_practice/metadata.tsv"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/bulk-comparison"))
    parser.add_argument("--species-map", type=Path, help="JSON mapping local compound labels to expected human/mouse/control")
    parser.add_argument("--exclude-samples", type=Path, help="Reviewed Sample 3/4 exclusions: CSV with sample_id and optional reason")
    args = parser.parse_args()
    species_map = {"231": "human", "3T3": "mouse", "human": "human", "mouse": "mouse", "water": "control", "blank": "control"}
    if args.species_map:
        species_map = json.loads(args.species_map.read_text())
    exclusions = None
    if args.exclude_samples:
        exclusions = pd.read_csv(args.exclude_samples, dtype={"sample_id": str})
        if ("sample_id" not in exclusions or exclusions.sample_id.isna().any()
                or exclusions.sample_id.str.strip().eq("").any() or exclusions.sample_id.duplicated().any()):
            parser.error("Exclusion CSV requires unique, non-missing sample_id values")
    compare(args.data_dir, args.bulk_expression, args.bulk_metadata, args.output_dir, species_map, exclusions)


if __name__ == "__main__":
    main()
