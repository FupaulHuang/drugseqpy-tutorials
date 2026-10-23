#!/usr/bin/env python3
"""Generate deterministic, entirely synthetic practice matrices (no source data)."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmwrite


def compressed(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as output:
        with gzip.GzipFile(filename="", mode="wb", fileobj=output, mtime=0) as handle:
            handle.write(payload)


def generate(root):
    seed = 271828
    rng = np.random.default_rng(seed)
    human_markers = ['VIM', 'FN1', 'CDH2', 'CD44', 'EGFR', 'MMP9', 'TGFB1', 'SNAI1', 'TWIST1', 'ITGB1', 'EPCAM', 'MT-ND1', 'MT-CO1', 'RPL3', 'RPS6']
    mouse_markers = ['Vim', 'Fn1', 'Col1a1', 'Col1a2', 'Acta2', 'mt-Nd1', 'mt-Co1', 'Rpl3', 'Rps6']
    summaries = []
    for number in range(1, 5):
        folder = root / f"sample_{number}"
        # Refuse to replace a user's own input bundle.
        if folder.exists():
            raise FileExistsError(f"{folder} already exists; choose a new --output directory")
        folder.mkdir(parents=True)
        if number < 3:
            genes = [f"ENSG{i:011d}" for i in range(1, 801)]
            symbols = human_markers + [f"DEMO_HUMAN_{i:04d}" for i in range(800 - len(human_markers))]
            groups = ['DMSO', 'Drug1', 'Drug2', 'water'] * 12
        else:
            genes = [f"ENSG{i:011d}" for i in range(1, 401)] + [f"ENSMUSG{i:011d}" for i in range(1, 401)]
            symbols = human_markers + [f"DEMO_HUMAN_{i:04d}" for i in range(400-len(human_markers))]
            symbols += mouse_markers + [f"DEMO_MOUSE_{i:04d}" for i in range(400-len(mouse_markers))]
            groups = ['231', '3T3'] * 23 + ['water', 'water']
        barcodes = [f"DEMO{i:04d}" for i in range(1, len(groups)+1)]
        # A fixed broad expression profile, mild plate shift, and well variation.
        baseline = np.random.default_rng(1729).lognormal(2.6, 0.8, len(genes))
        columns, metadata = [], []
        for well, (barcode, group) in enumerate(zip(barcodes, groups)):
            means = baseline.copy() * (1 + 0.06 * (number % 2))
            if number < 3:
                if group == 'Drug1':
                    means[30:100] *= 4
                    means[100:160] *= .25
                elif group == 'Drug2':
                    means[180:240] *= 2
                    means[240:280] *= .5
                elif group == 'water':
                    means *= .003
            else:
                if group == '231':
                    means[400:] *= .015
                elif group == '3T3':
                    means[:400] *= .02
                else:
                    means *= .002
                # Include a visibly mixed well and a low-depth well for review.
                if well == 8:
                    means = baseline.copy() * .45
                if well == 17:
                    means *= .004
            means *= rng.lognormal(0, .10)
            means *= rng.gamma(35, 1/35, len(genes))
            columns.append(rng.poisson(means).astype(np.int32))
            row, col = well % 8 + 1, well // 8 + 1
            metadata.append({
                'plate_id': f'Sample {number}', 'barcode': barcode,
                'well_id': f'Well-dT-{well+1}', 'compound': group,
                'dose': '0' if group in ('DMSO', 'water', '231', '3T3') else '1',
                'time': '24h', 'sample_type': 'media' if group == 'water' else 'DMSO' if group == 'DMSO' else 'treatment',
                'total_well': 96, 'x': col, 'y': row,
            })
        matrix = sparse.coo_matrix(np.column_stack(columns))
        buffer = io.BytesIO()
        mmwrite(buffer, matrix, field='integer', symmetry='general', comment='Entirely synthetic tutorial practice counts')
        compressed(folder/'matrix/matrix.mtx.gz', buffer.getvalue())
        compressed(folder/'matrix/features.tsv.gz', ''.join(f'{g}\t{s}\tGene Expression\n' for g,s in zip(genes,symbols)).encode())
        compressed(folder/'matrix/barcodes.tsv.gz', ''.join(f'{b}-1\n' for b in barcodes).encode())
        with (folder/'metadata.csv').open('w', newline='') as handle:
            writer=csv.DictWriter(handle, fieldnames=list(metadata[0]))
            writer.writeheader(); writer.writerows(metadata)
        summaries.append({'dataset':f'Sample {number}', 'wells':len(groups), 'features':len(genes), 'total_umi':int(matrix.sum())})
    checksums = {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(root.glob('sample_*/*')) if p.is_file()}
    checksums.update({str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(root.glob('sample_*/matrix/*'))})
    (root/'practice-manifest.json').write_text(json.dumps({'kind':'synthetic', 'seed':seed, 'datasets':summaries, 'sha256':checksums}, indent=2)+'\n')
    print(f'Generated four synthetic practice bundles in {root}')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('data'))
    args=parser.parse_args()
    generate(args.output)
