#!/usr/bin/env python3
"""Check the public API needed by the tutorial before running expensive steps."""
import inspect
import sys
from importlib.metadata import version

if sys.version_info < (3, 10):
    raise SystemExit('Use Python 3.10 or later.')
try:
    import drugseqpy as ds
except ImportError as error:
    raise SystemExit(f'Install drugseqpy and its dependencies first: {error}')

contracts = {
    'create_drugseq_from_files': ('metadata_key', 'min_overlap'),
    'merge_drugseq_objects': ('batch_key',),
    'compute_qc_metrics': ('mito_prefixes',),
    'build_plate_layout': ('nrow', 'ncol'),
    'generate_qc_report': ('star_root',),
    'assign_well_species': ('mapping',),
    'compute_mix_species_qc': ('assignment_threshold', 'min_species_umi'),
    'generate_mix_species_report': ('species_a_prefix', 'species_b_prefix'),
    'correct_batch': ('method', 'covariate_cols'),
    'compute_batch_effect_metrics': ('before', 'after'),
    'compute_multi_de': ('batch_cols',),
    'compute_cross_plate_reproducibility': ('reference',),
    'run_comparison_workflow': ('within_plate', 'gene_set_paths'),
    'run_enrichment': ('allow_network',),
}
missing = []
for name, parameters in contracts.items():
    function = getattr(ds, name, None)
    if not callable(function):
        missing.append(name)
    else:
        signature = inspect.signature(function)
        missing.extend(f'{name}({parameter}=...)' for parameter in parameters if parameter not in signature.parameters)
if missing:
    raise SystemExit('Install a compatible drugseqpy development revision. Missing API: ' + ', '.join(missing))
print(f'Python {sys.version.split()[0]}; drugseqpy {version("drugseqpy")}; tutorial API checks passed.')
print('API availability is checked; numerical reproducibility also requires recording the exact source revision and dependencies.')
