#!/usr/bin/env python3
"""Generate small synthetic bulk expression references without reading source data."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def generate(output: Path):
    if output.exists():
        raise FileExistsError(f"{output} already exists; choose a new --output directory")
    output.mkdir(parents=True)
    seed = 161803
    rng = np.random.default_rng(seed)
    genes = [f"ENSG{i:011d}" for i in range(1, 401)] + [f"ENSMUSG{i:011d}" for i in range(1, 401)]
    expression = pd.DataFrame(index=pd.Index(genes, name="feature_id"))
    metadata = []
    for species, compound, start in [("human", "231", 0), ("mouse", "3T3", 400)]:
        # Generated independently from the local practice count matrices.
        baseline = rng.lognormal(2.6, .9, 400)
        for study_number, unit in [(1, "counts"), (2, "FPKM")]:
            for replicate in (1, 2):
                ref = f"practice_{species}_{study_number}_{replicate}"
                values = baseline * rng.lognormal(0, .18, 400)
                values = rng.poisson(values * 20) if unit == "counts" else np.round(values, 6)
                expression[ref] = np.nan
                expression.loc[genes[start:start+400], ref] = values
                metadata.append({
                    "reference_id": ref, "study_id": f"practice_{species}_study_{study_number}",
                    "species": species, "compound": compound, "unit": unit, "primary": True,
                })
    expression.to_csv(output / "expression.tsv", sep="\t", na_rep="")
    pd.DataFrame(metadata).to_csv(output / "metadata.tsv", sep="\t", index=False)
    (output / "README.md").write_text(
        "# Synthetic bulk RNA-seq practice reference\n\n"
        "These eight references were generated independently from random distributions; "
        "they contain no observed biological data. The Ensembl-shaped identifiers are "
        "synthetic matching keys. The study labels describe simulation groups only.\n\n"
        "Human and mouse expression occupy separate rows. Blank cells mean unmeasured, "
        "not zero. Counts and FPKM columns retain distinct units. "
        "Use the dataset to exercise file formats and correlation calculations, "
        "not to draw conclusions about biology or methods.\n"
    )
    checksums = {name: hashlib.sha256((output / name).read_bytes()).hexdigest()
                 for name in ("expression.tsv", "metadata.tsv", "README.md")}
    (output / "manifest.json").write_text(json.dumps({
        "kind": "entirely synthetic", "seed": seed,
        "features": 800, "references": 8, "sha256": checksums,
    }, indent=2) + "\n")
    print(f"Generated synthetic bulk practice references in {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("data/bulk_practice"))
    generate(parser.parse_args().output)
