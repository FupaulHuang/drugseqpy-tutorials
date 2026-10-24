#!/usr/bin/env python3
"""Prepare an explicitly selected, anonymized Drug-seq count subset.

Inputs stay outside the tutorial repository. Metadata CSV (optionally .gz)
requires plate_barcode, well_index, plate_well, compound, treatment_dose,
treatment_dose_unit and timepoint. A manifest CSV supplies label (Sample 1,
Sample 2, ...), source_plate and counts (a genes-by-wells TSV, optionally .gz).
Counts columns must be source_plate + '_' + well_index; its first column
contains unique feature IDs. A gene-symbol/contig composite ID can be retained;
gene_symbol is the text before its final comma.

Selection JSON example:
{"compounds": {"Treatment A": "Drug1", "Treatment B": "Drug2"},
 "doses": [0.1, 1.0], "dose_unit": "uM", "time": 24,
 "control": "DMSO"}

Select plates and conditions before examining expression outcomes. All control
wells at the specified time are retained. This script does not exclude wells
on QC grounds, normalize counts, collapse gene symbols or fill missing genes.
The control setting identifies the source label for DMSO vehicle wells and
exports it as DMSO; it must not be used to relabel a different vehicle.
If cell.line.name is supplied, all selected wells must share one complete value.
Source identifiers and locations appear only in an optional private audit.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


def compressed(path: Path, payload: bytes):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as target:
        with gzip.GzipFile(filename="", fileobj=target, mode="wb", mtime=0) as stream:
            stream.write(payload)


def prepare(metadata_path: Path, manifest_path: Path, selection_path: Path,
            output: Path, private_audit: Path | None = None):
    if output.exists():
        raise FileExistsError(f"{output} exists; choose a new output directory")
    metadata = pd.read_csv(metadata_path, dtype=str, keep_default_na=False)
    required = ["plate_barcode", "well_index", "plate_well", "compound",
                "treatment_dose", "treatment_dose_unit", "timepoint"]
    if set(required) - set(metadata):
        raise ValueError(f"Metadata needs columns: {sorted(set(required) - set(metadata))}")
    if metadata[required].apply(lambda s: s.str.strip().eq("").any()).any():
        raise ValueError("Required metadata fields must be complete")
    manifest = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    if set(["label", "source_plate", "counts"]) - set(manifest) or manifest.empty:
        raise ValueError("Manifest needs at least one row with label, source_plate and counts")
    if manifest[["label", "source_plate", "counts"]].eq("").any().any():
        raise ValueError("Manifest entries must be complete")
    if not manifest.label.is_unique or not manifest.source_plate.is_unique:
        raise ValueError("Manifest labels and source plates must each be unique")
    if not manifest.label.str.fullmatch(r"Sample [1-9]\d*").all():
        raise ValueError("Use neutral manifest labels Sample 1, Sample 2, ...")
    selection = json.loads(selection_path.read_text())
    compounds = selection.get("compounds", {})
    if not isinstance(compounds, dict) or not compounds or any(not str(k).strip() for k in compounds):
        raise ValueError("Selection compounds must map source names to neutral aliases")
    aliases = list(compounds.values())
    if len(set(aliases)) != len(aliases) or not all(re.fullmatch(r"Drug[1-9]\d*", str(a)) for a in aliases):
        raise ValueError("Compound aliases must be distinct Drug1, Drug2, ... labels")
    doses = np.asarray(selection.get("doses", []), dtype=float)
    if not len(doses) or not np.isfinite(doses).all() or (doses < 0).any() or len(set(doses)) != len(doses):
        raise ValueError("Selection doses must be distinct finite nonnegative values")
    unit = str(selection.get("dose_unit", ""))
    time = float(selection["time"])
    if not unit or not np.isfinite(time) or time < 0:
        raise ValueError("Selection needs a dose_unit and a finite nonnegative time")
    control = str(selection.get("control", "DMSO"))
    if control in compounds:
        raise ValueError("Select the control separately from treatment compounds")
    dose_values = pd.to_numeric(metadata.treatment_dose, errors="raise")
    time_values = pd.to_numeric(metadata.timepoint, errors="raise")
    if not np.isfinite(dose_values).all() or not np.isfinite(time_values).all():
        raise ValueError("Dose and time metadata must be finite numbers")
    metadata = metadata.assign(_dose=dose_values, _time=time_values)
    loaded, audits, designs, cell_lines = [], [], [], []
    expected_conditions = {"DMSO"} | {f"{alias}_{dose:g}{unit}" for alias in aliases for dose in doses}
    for entry in manifest.itertuples(index=False):
        plate_meta = metadata.loc[metadata.plate_barcode.eq(entry.source_plate)].copy()
        if plate_meta.empty:
            raise ValueError(f"No metadata found for manifest label {entry.label}")
        keys = plate_meta.plate_barcode + "_" + plate_meta.well_index
        if not keys.is_unique or not plate_meta.plate_well.is_unique:
            raise ValueError(f"{entry.label}: duplicate barcode or physical well")
        selected = plate_meta.loc[
            plate_meta._time.eq(time) & (
                plate_meta.compound.eq(control)
                | (plate_meta.compound.isin(compounds) & plate_meta._dose.isin(doses))
            )
        ].copy()
        if "cell.line.name" in selected:
            cells = selected["cell.line.name"]
            if cells.isna().any() or cells.str.strip().eq("").any():
                raise ValueError(f"{entry.label}: selected cell.line.name values must be complete")
            cell_lines.extend(cells.to_list())
        if not selected.treatment_dose_unit.eq(unit).all():
            raise ValueError(f"{entry.label}: selected dose units disagree with the selection")
        if not selected.loc[selected.compound.eq(control), "_dose"].eq(0).all():
            raise ValueError(f"{entry.label}: expected a recorded zero-dose control")
        positions = selected.plate_well.str.upper().str.extract(r"^([A-P])(\d{1,2})$")
        if positions.isna().any().any():
            raise ValueError(f"{entry.label}: use unambiguous A01–P24 physical well labels")
        columns = pd.to_numeric(positions[1])
        if not columns.between(1, 24).all():
            raise ValueError(f"{entry.label}: physical well column is outside 1–24")
        selected["_row"] = positions[0].map(lambda s: ord(s) - ord("A") + 1)
        selected["_column"] = columns
        selected = selected.sort_values(["_row", "_column"])
        selected["_key"] = selected.plate_barcode + "_" + selected.well_index
        selected["_well"] = selected.apply(lambda row: f"{chr(int(row['_row']) + 64)}{int(row['_column']):02d}", axis=1)
        if not selected._well.is_unique:
            raise ValueError(f"{entry.label}: physical well aliases are duplicated")
        public = pd.DataFrame({
            "plate_id": entry.label, "barcode": selected._well,
            "sample_id": entry.label + "_" + selected._well,
            "well_id": selected._well,
            "compound": selected.compound.map(compounds).fillna("DMSO"),
            "dose": selected._dose, "dose_unit": unit, "time": f"{time:g}h",
            "sample_type": np.where(selected.compound.eq(control), "DMSO", "treatment"),
            "total_well": 384, "x": selected._column, "y": selected._row,
        })
        public["condition"] = ["DMSO" if drug == "DMSO" else f"{drug}_{dose:g}{unit}"
                               for drug, dose in zip(public.compound, public.dose)]
        design = public.groupby("condition").size()
        if set(design.index) != expected_conditions or (design < 2).any():
            raise ValueError(f"{entry.label}: every requested treatment-dose and DMSO needs >=2 wells")
        designs.append(design)
        path = Path(entry.counts)
        if not path.is_absolute():
            path = manifest_path.parent / path
        header = pd.read_csv(path, sep="\t", nrows=0).columns
        if len(header) < 2 or not header[1:].is_unique:
            raise ValueError(f"{entry.label}: invalid count-table header")
        if set(header[1:]) != set(keys):
            raise ValueError(f"{entry.label}: count columns and full plate metadata keys must agree exactly")
        counts = pd.read_csv(path, sep="\t", usecols=[header[0]] + selected._key.tolist(), index_col=0)
        if counts.index.isna().any():
            raise ValueError(f"{entry.label}: feature IDs must be present")
        counts.index = counts.index.astype(str)
        if not counts.index.is_unique or (counts.index.str.strip() == "").any():
            raise ValueError(f"{entry.label}: feature IDs must be present and unique")
        counts = counts.loc[:, selected._key]
        values = counts.to_numpy()
        if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
            raise ValueError(f"{entry.label}: provide finite nonnegative integer UMI counts")
        if (values.sum(axis=0) == 0).any():
            raise ValueError(f"{entry.label}: a selected well has no observed UMIs")
        counts.columns = public.barcode.to_list()
        loaded.append((entry.label, counts, public.reset_index(drop=True)))
        audits.append({
            "label": entry.label, "source_plate": entry.source_plate,
            "count_path": str(path.resolve()), "input_features": len(counts),
            "full_plate_wells": len(keys), "selected_wells": len(selected),
            "exact_full_metadata_count_key_set": True,
            "metadata_original_order_matches_counts": list(keys) == list(header[1:]),
            "well_map": selected[["_key", "_well", "compound", "_dose"]].to_dict(orient="records"),
        })
    if not all(designs[0].sort_index().equals(d.sort_index()) for d in designs[1:]):
        raise ValueError("Selected plates must have equal counts for each requested condition")
    if cell_lines and len(set(cell_lines)) != 1:
        raise ValueError("All selected plates and conditions must use the same cell.line.name")
    shared = set.intersection(*(set(counts.index) for _, counts, _ in loaded))
    features = pd.Index(sorted(shared))
    if len(features) < 3:
        raise ValueError("Fewer than three shared features")
    output.mkdir(parents=True)
    public_summary = []
    for (label, counts, public), audit in zip(loaded, audits):
        missing = counts.index.difference(features)
        dropped_umi = int(counts.loc[missing].to_numpy().sum())
        audit.update({"dropped_features": missing.to_list(), "dropped_umi_selected_wells": dropped_umi})
        counts = counts.loc[features]
        folder = output / ("sample_" + label.split()[-1])
        matrix = sparse.coo_matrix(counts.to_numpy(dtype=np.int64))
        buffer = io.BytesIO()
        mmwrite(buffer, matrix, field="integer", symmetry="general",
                comment="Anonymized selected Drug-seq UMI counts; features x wells")
        compressed(folder / "matrix/matrix.mtx.gz", buffer.getvalue())
        symbols = features.str.rsplit(",", n=1).str[0]
        feature_text = "".join(
            f"{feature}\t{symbol}\t{'Spike-in' if feature.startswith('ERCC-') else 'Gene Expression'}\n"
            for feature, symbol in zip(features, symbols)
        )
        compressed(folder / "matrix/features.tsv.gz", feature_text.encode())
        compressed(folder / "matrix/barcodes.tsv.gz", "".join(f"{b}\n" for b in public.barcode).encode())
        public.to_csv(folder / "metadata.csv", index=False)
        qc = public[["sample_id", "plate_id", "well_id", "compound", "dose", "condition"]].copy()
        qc["total_umi"] = np.asarray(matrix.sum(axis=0)).ravel()
        qc["n_genes_detected"] = np.asarray((matrix > 0).sum(axis=0)).ravel()
        qc.to_csv(folder / "input_depth_audit.csv", index=False)
        public_summary.append({
            "plate_id": label, "selected_wells": len(public), "shared_features": len(features),
            "dropped_nonshared_features": len(missing), "dropped_umi_selected_wells": dropped_umi,
            "condition_counts": public.condition.value_counts().sort_index().to_dict(),
        })
    checksums = {str(path.relative_to(output)): hashlib.sha256(path.read_bytes()).hexdigest()
                 for path in sorted(output.glob("sample_*/*/*")) if path.is_file()}
    checksums.update({str(path.relative_to(output)): hashlib.sha256(path.read_bytes()).hexdigest()
                      for path in sorted(output.glob("sample_*/*.csv"))})
    (output / "preparation_summary.json").write_text(json.dumps({
        "kind": "selected experimental UMI counts", "plates": public_summary,
        "feature_policy": "shared feature IDs explicitly sorted; no zero filling or symbol collapsing",
        "selection_policy": "prespecified plates, compounds and doses; all matching controls; no QC exclusions",
        "sha256": checksums,
    }, indent=2) + "\n")
    if private_audit:
        private_audit.parent.mkdir(parents=True, exist_ok=True)
        private_audit.write_text(json.dumps({"selection": selection, "plates": audits}, indent=2) + "\n")
    print(f"Prepared {len(loaded)} plates, {sum(len(obs) for _, _, obs in loaded)} wells, {len(features)} shared features in {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--private-audit", type=Path, help="Optional private source-to-alias audit; keep outside a public repository")
    args = parser.parse_args()
    prepare(args.metadata, args.manifest, args.selection, args.output_dir, args.private_audit)


if __name__ == "__main__":
    main()
