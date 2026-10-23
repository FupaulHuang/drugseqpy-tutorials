#!/usr/bin/env python3
"""Compare expression profiles among wells on one plate.

Heatmaps use genes detected anywhere in their displayed well set. Pair
statistics use genes detected in either member of that pair. All statistics
use exact log10(raw UMI + 1); optional display jitter never enters a statistic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata

import drugseqpy as ds


def natural_key(value):
    return tuple(int(part) if part.isdigit() else part.casefold()
                 for part in re.split(r"(\d+)", str(value)))


def correlation(x, y):
    """Return Pearson correlation, without assigning values to constant data."""
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    x = x - x.mean()
    y = y - y.mean()
    denominator = np.linalg.norm(x) * np.linalg.norm(y)
    if denominator == 0:
        return np.nan
    return float(np.clip(np.dot(x, y) / denominator, -1.0, 1.0))


def pair_statistics(first, second):
    """Drop jointly undetected genes, then calculate both correlation measures."""
    keep = (first > 0) | (second > 0)
    x, y = np.log10(first[keep] + 1), np.log10(second[keep] + 1)
    pearson = correlation(x, y)
    spearman = correlation(rankdata(x), rankdata(y))
    status = "ok" if np.isfinite(pearson) and np.isfinite(spearman) else "undefined"
    if status == "ok":
        slope, intercept = np.polyfit(x, y, 1)
    else:
        slope, intercept = np.nan, np.nan
    return {
        "n_genes_pair": int(keep.sum()), "pearson_r": pearson,
        "spearman_rho": spearman, "fit_slope": float(slope),
        "fit_intercept": float(intercept), "status": status,
    }, x, y


def matrix_correlations(counts, method):
    """Use one detected-gene universe for the entire displayed heatmap."""
    keep = (counts > 0).any(axis=0)
    values = np.log10(counts[:, keep] + 1)
    if method == "spearman":
        values = np.asarray([rankdata(row) for row in values])
    result = np.full((len(counts), len(counts)), np.nan)
    for i in range(len(counts)):
        for j in range(i, len(counts)):
            result[i, j] = result[j, i] = correlation(values[i], values[j])
    return result, int(keep.sum())


def save_figure(fig, path):
    fig.savefig(path, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def draw_heatmap(counts, obs, conditions, flagged, method, title, path):
    result, n_genes = matrix_correlations(counts, method)
    n = len(obs)
    fig, ax = plt.subplots(figsize=(10.5, 9.5) if n > 16 else (8.4, 7.5))
    cmap = plt.get_cmap("RdYlBu_r").copy()
    cmap.set_bad("#d0d4d7")
    im = ax.imshow(np.ma.masked_invalid(result), cmap=cmap, vmin=-1, vmax=1)
    labels = obs["well_id"].astype(str).tolist()
    size = 5.5 if n > 32 else 7.5
    ax.set_xticks(np.arange(n), labels, rotation=90, fontsize=size)
    ax.set_yticks(np.arange(n), labels, fontsize=size)
    for ticks in (ax.get_xticklabels(), ax.get_yticklabels()):
        for tick, sample_id in zip(ticks, obs["sample_id"]):
            if sample_id in flagged:
                tick.set_color("#b4233b")
                tick.set_fontweight("bold")
    # Group borders preserve the condition ordering; no clustering is implied.
    group_keys = list(obs[conditions].itertuples(index=False, name=None))
    starts = [0] + [i for i in range(1, n) if group_keys[i] != group_keys[i - 1]]
    ends = starts[1:] + [n]
    for start, end in zip(starts, ends):
        if start:
            ax.axhline(start - 0.5, color="white", lw=1.4)
            ax.axvline(start - 0.5, color="white", lw=1.4)
        if len(starts) <= 6:
            label = str(obs.iloc[start]["compound"])
            ax.text((start + end - 1) / 2, 1.012, label,
                    transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                    fontsize=8, fontweight="bold")
    ax.set_xlabel("Wells ordered by condition, then well number", labelpad=13)
    ax.set_ylabel("Wells in the same order")
    label = "Pearson r" if method == "pearson" else "Spearman rho"
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.035, label=label)
    ax.set_title(f"{title}\n{n:,} wells · {n_genes:,} genes detected in this well set",
                 pad=17, fontsize=12, fontweight="bold")
    fig.text(0.5, 0.015,
             "log10(raw UMI + 1) · grey = undefined · red labels = flagged, retained",
             ha="center", fontsize=8, color="#536574")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save_figure(fig, path)
    return result, n_genes


def draw_scatter(row, x, y, path, synthetic, jitter):
    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    if jitter:
        seed = int(hashlib.sha256(row["pair_id"].encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(seed)
        dx = np.maximum(0, x + rng.uniform(-jitter, jitter, len(x)))
        dy = np.maximum(0, y + rng.uniform(-jitter, jitter, len(y)))
    else:
        dx, dy = x, y
    limit = max(1.0, float(np.ceil(max(x.max(), y.max()))))
    ax.scatter(dx, dy, s=5, alpha=0.42, color="#167d91", edgecolors="none", rasterized=True)
    ax.plot([0, limit], [0, limit], "--", color="#6a7884", lw=1.1, label="Identity: y = x")
    ax.plot([0, limit], [row["fit_intercept"], row["fit_intercept"] + row["fit_slope"] * limit],
            color="#d55e00", lw=1.3, label="OLS fit: y on x")
    ax.set(xlim=(0, limit), ylim=(0, limit), aspect="equal",
           xlabel=f'{row["well_a"]}: log10(raw UMI + 1)',
           ylabel=f'{row["well_b"]}: log10(raw UMI + 1)')
    context = f'{row["plate_id"]} · {row["condition_label"]}'
    if synthetic:
        context = "Synthetic practice · " + context
    ax.set_title(f'{context}\n{row["well_a"]} versus {row["well_b"]}',
                 fontsize=11, fontweight="bold")
    note = (f'Pearson r = {row["pearson_r"]:.4f}\n'
            f'Spearman rho = {row["spearman_rho"]:.4f}\n'
            f'{row["n_genes_pair"]:,} genes detected in either well\n'
            f'QC flags: {row["flagged_a"]} / {row["flagged_b"]}')
    ax.text(0.025, 0.975, note, transform=ax.transAxes, ha="left", va="top", fontsize=8.5,
            bbox={"facecolor": "white", "edgecolor": "#b6c1c8", "alpha": 0.94, "pad": 5})
    ax.legend(loc="lower right", frameon=False, fontsize=8)
    ax.grid(color="#e3e9ec", lw=0.65)
    ax.set_axisbelow(True)
    if jitter:
        fig.text(0.5, 0.015, f"Display jitter ±{jitter:g}; statistics use unjittered values",
                 ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.025, 1, 1))
    save_figure(fig, path)


def choose_pairs(pairs, explicit, cap, control):
    valid = pairs[pairs["status"].eq("ok")].copy()
    if explicit:
        missing = set(explicit) - set(valid["pair_id"])
        if missing:
            raise ValueError(f"Requested pair IDs are absent or undefined: {sorted(missing)}")
        if len(explicit) > cap:
            raise ValueError("Increase --max-scatters to include every requested --pair.")
        return [(pair_id, "explicit request") for pair_id in dict.fromkeys(explicit)]
    if cap == 0 or valid.empty:
        return []
    selected = []
    def add(frame, reason, descending=False):
        frame = frame[~frame["pair_id"].isin([item[0] for item in selected])]
        if not frame.empty and len(selected) < cap:
            item = frame.sort_values(["pearson_r", "pair_id"], ascending=[not descending, True]).iloc[0]
            selected.append((item["pair_id"], reason))
    # Deliberately illustrate the spread; these examples are not a representative sample.
    add(valid[valid["compound"].eq(control)], "lowest Pearson control pair")
    treatment = valid[~valid["compound"].eq(control)
                      & ~valid["compound"].str.casefold().isin(["water", "blank", "media"])]
    add(treatment, "lowest Pearson treatment pair")
    add(treatment, "highest Pearson treatment pair", descending=True)
    while len(selected) < min(cap, len(valid)):
        add(valid, "remaining pairs in ascending Pearson order")
    return selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--sample", default="sample_1", help="One sample subdirectory containing matrix/ and metadata.csv")
    parser.add_argument("--output-dir", type=Path, default=Path("results/well-repeatability"))
    parser.add_argument("--condition-columns", nargs="+", default=["compound", "dose", "time"])
    parser.add_argument("--heatmap-method", choices=["pearson", "spearman"], default="pearson")
    parser.add_argument("--max-scatters", type=int, default=6, help="Cap figures only; all pairs are still computed")
    parser.add_argument("--pair", action="append", default=[], help="Explicit pair_id from replicate_pairs.csv; repeat to select several")
    parser.add_argument("--control", default="DMSO", help="Control label used only for automatic illustration selection")
    parser.add_argument("--flagged-samples", type=Path, help="CSV with sample_id column; mark these wells without removing them")
    parser.add_argument("--display-jitter", type=float, default=0.0, help="Optional uniform scatter jitter in log10 units; statistics unchanged")
    parser.add_argument("--synthetic", action="store_true", help="Label figures as synthetic practice illustrations")
    args = parser.parse_args()
    if args.max_scatters < 0 or not 0 <= args.display_jitter <= 0.05:
        parser.error("--max-scatters must be nonnegative; --display-jitter must be between 0 and 0.05")
    condition_columns = list(dict.fromkeys(args.condition_columns))
    if "compound" not in condition_columns or "plate_id" in condition_columns:
        parser.error("Include compound in --condition-columns; plate_id is always checked separately")
    source = args.data_dir / args.sample
    obj = ds.create_drugseq_from_files(source / "matrix", source / "metadata.csv",
                                      metadata_key="barcode", min_overlap=1.0)
    obj.adata.obs_names = obj.obs["sample_id"].astype(str)
    if not obj.adata.var_names.is_unique:
        raise ValueError("Feature identifiers must be unique before profile comparison")
    obs = obj.obs.copy().reset_index(drop=True)
    required = ["plate_id", "sample_id", "well_id"] + condition_columns
    missing = set(required) - set(obs)
    if missing:
        raise ValueError(f"Missing metadata columns: {sorted(missing)}")
    if obs[required].isna().any().any():
        raise ValueError("Condition fields and well identifiers must be complete")
    for column in required:
        obs[column] = obs[column].astype(str)
        if obs[column].str.strip().eq("").any():
            raise ValueError(f"Empty values in {column}")
    if obs["plate_id"].nunique() != 1 or len(obs) < 2:
        raise ValueError("This script needs at least two wells from exactly one plate")
    if not obs["sample_id"].is_unique or not obs["well_id"].is_unique:
        raise ValueError("sample_id and well_id must uniquely identify wells on this plate")
    counts = obj.layers["counts"]
    if counts.shape[0] * counts.shape[1] * 8 > 512 * 1024**2:
        raise ValueError("Dense count matrix would exceed 512 MiB; process a smaller well set")
    counts = counts.toarray() if sparse.issparse(counts) else np.asarray(counts)
    counts = counts.astype(float)
    if not np.isfinite(counts).all() or (counts < 0).any() or not np.equal(counts, np.floor(counts)).all():
        raise ValueError("The counts layer must contain finite, nonnegative integer UMIs")
    order = sorted(range(len(obs)), key=lambda i: tuple(obs.iloc[i][condition_columns])
                   + (natural_key(obs.iloc[i]["well_id"]),))
    obs, counts = obs.iloc[order].reset_index(drop=True), counts[order]
    flagged = set()
    if args.flagged_samples:
        flags = pd.read_csv(args.flagged_samples, dtype=str)
        if "sample_id" not in flags or flags["sample_id"].isna().any():
            raise ValueError("Flag table requires a nonmissing sample_id column")
        flagged = set(flags["sample_id"])
        if flagged - set(obs["sample_id"]):
            raise ValueError("Flag table contains sample IDs absent from the input plate")
    output = args.output_dir
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    obs["flagged_for_review"] = obs["sample_id"].isin(flagged)
    obs.to_csv(output / "well_order.csv", index=False)
    prefix = "Synthetic practice · " if args.synthetic else ""
    plate = obs["plate_id"].iloc[0]
    heatmap, heatmap_n_genes = draw_heatmap(
        counts, obs, condition_columns, flagged, args.heatmap_method,
        f"{prefix}{plate} · all-well profile similarity", figures / "all_wells_heatmap.png")
    pd.DataFrame(heatmap, index=obs["sample_id"], columns=obs["sample_id"]).to_csv(output / "all_wells_correlations.csv")
    pair_rows, group_rows = [], []
    for group_number, (key, indexes) in enumerate(obs.groupby(condition_columns, sort=True).groups.items(), 1):
        if not isinstance(key, tuple):
            key = (key,)
        group = dict(zip(condition_columns, key))
        indexes = list(indexes)
        group_id = f"condition-{group_number:03d}"
        label = "; ".join(f"{name}={value}" for name, value in group.items())
        group_rows.append({"condition_id": group_id, **group, "n_wells": len(indexes), "condition_label": label})
        if len(indexes) >= 2:
            corr, n_genes = draw_heatmap(
                counts[indexes], obs.iloc[indexes], condition_columns, flagged, args.heatmap_method,
                f"{prefix}{plate}\n{label}", figures / f"{group_id}_heatmap.png")
            pd.DataFrame(corr, index=obs.iloc[indexes]["sample_id"], columns=obs.iloc[indexes]["sample_id"]).to_csv(output / f"{group_id}_correlations.csv")
            group_rows[-1]["heatmap_n_genes"] = n_genes
        for i, j in combinations(indexes, 2):
            stats, _, _ = pair_statistics(counts[i], counts[j])
            pair_rows.append({
                "pair_id": f"pair-{len(pair_rows) + 1:05d}", "plate_id": plate,
                "condition_id": group_id, **group, "condition_label": label,
                "sample_a": obs.iloc[i]["sample_id"], "sample_b": obs.iloc[j]["sample_id"],
                "well_a": obs.iloc[i]["well_id"], "well_b": obs.iloc[j]["well_id"],
                "flagged_a": obs.iloc[i]["sample_id"] in flagged,
                "flagged_b": obs.iloc[j]["sample_id"] in flagged,
                "n_wells_condition": len(indexes), **stats,
            })
    if not pair_rows:
        raise ValueError("No condition has at least two wells; inspect your condition columns")
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(output / "replicate_pairs.csv", index=False)
    summary = pairs.groupby(["condition_id"] + condition_columns, dropna=False).agg(
        n_wells=("n_wells_condition", "first"), n_pairs=("pair_id", "size"),
        n_defined=("pearson_r", "count"), pearson_median=("pearson_r", "median"),
        pearson_min=("pearson_r", "min"), pearson_max=("pearson_r", "max"),
        spearman_median=("spearman_rho", "median"),
        n_genes_pair_min=("n_genes_pair", "min"), n_genes_pair_max=("n_genes_pair", "max"),
    ).reset_index()
    summary.to_csv(output / "condition_summary.csv", index=False)
    pd.DataFrame(group_rows).to_csv(output / "condition_index.csv", index=False)
    chosen = choose_pairs(pairs, args.pair, args.max_scatters, args.control)
    by_sample = {sample_id: index for index, sample_id in enumerate(obs["sample_id"])}
    selected_rows = []
    for number, (pair_id, reason) in enumerate(chosen, 1):
        row = pairs.set_index("pair_id", drop=False).loc[pair_id].to_dict()
        _, x, y = pair_statistics(counts[by_sample[row["sample_a"]]], counts[by_sample[row["sample_b"]]])
        name = f"scatter-{number:02d}.png"
        draw_scatter(row, x, y, figures / name, args.synthetic, args.display_jitter)
        selected_rows.append({"pair_id": pair_id, "selection_reason": reason, "figure": name})
    pd.DataFrame(selected_rows, columns=["pair_id", "selection_reason", "figure"]).to_csv(output / "selected_pairs.csv", index=False)
    settings = {
        "plate_id": plate, "synthetic_illustration": args.synthetic,
        "n_wells": len(obs), "n_input_genes": counts.shape[1], "all_wells_heatmap_n_genes": heatmap_n_genes,
        "n_pairs": len(pairs), "n_defined_pairs": int(pairs["status"].eq("ok").sum()),
        "n_flagged_wells_retained": len(flagged), "condition_columns": ["plate_id"] + condition_columns,
        "transformation": "log10(raw UMI + 1)", "heatmap_method": args.heatmap_method,
        "heatmap_gene_universe": "detected in at least one well of the displayed set",
        "pair_gene_universe": "detected in at least one member of the pair",
        "pair_statistics": ["Pearson r", "Spearman rho with average ranks"],
        "undefined_policy": "fewer than three retained genes or constant profile: NaN, no scatter",
        "display_jitter": args.display_jitter, "statistics_use_jitter": False,
        "n_scatter_figures": len(chosen), "max_scatters": args.max_scatters,
        "no_wells_excluded": True,
    }
    (output / "run_settings.json").write_text(json.dumps(settings, indent=2) + "\n")
    print(f"Saved {len(obs)}-well heatmap, condition heatmaps, {len(pairs)} pair rows, and {len(chosen)} scatters to {output}")


if __name__ == "__main__":
    main()
