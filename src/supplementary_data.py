"""Supplementary Data 1 and 2: Class I phosphosite matrices for the DVP datasets.

Supplementary Data 1 = FFPE lung adenocarcinoma, all microdissected samples (Fig. 4h-l).
Supplementary Data 2 = phosphoDVP mouse-brain neuron populations, main experiment (Fig. 5e-f).

Both are built from the wide Spectronaut PTM Site Reports with the same per-run Class I
processing as the figures (``core.process_ptm_site_report``, cutoff 0.75). All runs are
kept (no sample QC), so the unique-site totals match the numbers stated in the text.
One row per site and multiplicity (the quantitative entity); the unique-site total
collapses multiplicity, as everywhere else in the manuscript.

Run from the repo root:  python src/supplementary_data.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from core import process_ptm_site_report

CLASS_I_CUTOFF = 0.75
OUT_DIR = Path('tables')
PROC_META = {'Protein_group', 'Gene_group', 'PTM_0_aa', 'PTM_pos', 'PTM_mult123',
             'PTM_flank', 'PTM_Collapse_key', 'PTM_localization', 'UPD_seq'}
LUNG_REPORT = Path('pride_data/analysis_data/revision/figure4')
LUNG_LEGEND = Path('tables/Supplementary_Table_3.xlsx')
NEURON_REPORT = Path('pride_data/analysis_data/revision/figure5')
NEURON_COND = Path('data/pDVP_neurons_conditionSetup.tsv')
LUNG_TYPE = {'tumor (adenocarcinoma)': 'AC', 'normal lung': 'H'}


def lung_samples(runs):
    legend = pd.read_excel(LUNG_LEGEND).set_index('384_well_position')
    rows = []
    for run in runs:
        # one run carries a timestamp suffix after the well (..._phosphoDVP_K13_20250924081708)
        well = re.search(r'phosphoDVP_([A-P]\d+)', run).group(1)
        if well not in legend.index:
            # O15/O17/O19: acquired but not in the sample legend (23-56 sites each, not tissue samples)
            continue
        info = legend.loc[well]
        rows.append({'Sample': f"{LUNG_TYPE[info['sample_type']]}_{info['patient']}_{well}",
                     'Raw file': f'{run}.raw', 'Well (384)': well, 'Tissue': info['sample_type'],
                     'Patient': info['patient'], 'Area (x1000 µm²)': info['area (*1000, sq um)']})
    return pd.DataFrame(rows)


def neuron_samples(runs):
    cond = pd.read_csv(NEURON_COND, sep='\t')
    cond['run'] = cond['Run Label'].str.replace(r'\.raw$', '', regex=True)
    cond = cond.set_index('run')
    rows = []
    for run in runs:
        c, rep = cond.loc[run, 'Condition'], cond.loc[run, 'Replicate']
        population = re.sub(r'\d+$', '', c)
        rows.append({'Sample': f'{c}_rep{rep}', 'Raw file': f'{run}.raw', 'Population': population,
                     'Cell type': 'excitatory' if population.startswith('act') else 'inhibitory',
                     'Region': population.split('_')[1], 'Condition': c, 'Replicate': int(rep)})
    return pd.DataFrame(rows)


def site_matrix(report_path, annotate):
    sd = process_ptm_site_report(pd.read_csv(report_path, sep='\t', low_memory=False),
                                 cutoff=CLASS_I_CUTOFF)['site_data']
    samples = annotate([c for c in sd.columns if c not in PROC_META])
    samples = samples.sort_values('Sample', kind='stable').reset_index(drop=True)
    run_to_sample = dict(zip(samples['Raw file'].str.replace(r'\.raw$', '', regex=True), samples['Sample']))
    runs = list(run_to_sample)
    sd = sd[sd[runs].notna().any(axis=1)].reset_index(drop=True)
    gene = sd['Gene_group'].where(~sd['Gene_group'].isin(['', 'nan']), '')
    sites = pd.DataFrame({
        'Protein': sd['Protein_group'],
        'Gene': gene,
        'Site': (gene.where(gene != '', sd['Protein_group']) + '_' + sd['PTM_0_aa']
                 + sd['PTM_pos'].astype(str)),
        'Residue': sd['PTM_0_aa'],
        'Position': sd['PTM_pos'],
        'Multiplicity': sd['PTM_mult123'].map({1: '1', 2: '2', 3: '3+'}),
        'Flanking sequence': sd['PTM_flank'],
        'Best localization probability': sd['PTM_localization'].round(4),
        'Class I runs': sd[runs].notna().sum(axis=1),
    })
    quant = sd[runs].rename(columns=run_to_sample)[samples['Sample']]
    # gene-less entries ('') sort after named genes
    out = pd.concat([sites, quant], axis=1).sort_values(
        ['Gene', 'Protein', 'Position', 'Multiplicity'],
        key=lambda s: s.replace('', '~') if s.name == 'Gene' else s)
    return out.reset_index(drop=True), samples


def readme(title, source, sites, samples):
    n_sites = sites['Site'].nunique()
    return pd.DataFrame({'Supplementary Data': [
        title,
        f'Source: Spectronaut PTM Site Report {source.name} (PRIDE deposit).',
        f'{n_sites:,} unique Class I phosphosites (multiplicity collapsed) across {len(samples)} runs; '
        f'{len(sites):,} rows (one per site and multiplicity).',
        f'Class I: localization probability >= {CLASS_I_CUTOFF} in the given run; a run value is reported '
        'only where the site is Class I in that run, otherwise blank.',
        'Quantification: log2 Spectronaut PTM.Quantity (site-determining fragments), as exported, without '
        'parent-protein normalization or imputation.',
        'Best localization probability: maximum across runs. Class I runs: number of runs with a value.',
        'Sheet "Samples" maps the column names to raw files and sample annotation.',
    ]})


def write(path, title, source, sites, samples):
    with pd.ExcelWriter(path) as xl:
        readme(title, source, sites, samples).to_excel(xl, sheet_name='README', index=False)
        sites.to_excel(xl, sheet_name='Class I sites', index=False)
        samples.to_excel(xl, sheet_name='Samples', index=False)
    print(f'{path}: {sites["Site"].nunique():,} unique Class I sites, {len(sites):,} rows, '
          f'{len(samples)} runs')


def main():
    lung_src = next(LUNG_REPORT.glob('*nanoPhos_lungAC_Report.tsv'))
    sites, samples = site_matrix(lung_src, lung_samples)
    write(OUT_DIR / 'Supplementary_Data_1.xlsx',
          'Supplementary Data 1. Class I phosphosites in laser-microdissected FFPE lung adenocarcinoma (AC) '
          'and adjacent healthy-appearing lung (H), all microdissected samples (Fig. 4h-l).',
          lung_src, sites, samples)

    neuron_src = next(NEURON_REPORT.glob('*Whisper80_Report.tsv'))
    sites, samples = site_matrix(neuron_src, neuron_samples)
    write(OUT_DIR / 'Supplementary_Data_2.xlsx',
          'Supplementary Data 2. Class I phosphosites in laser-microdissected excitatory and inhibitory '
          'neurons from mouse cortex and subcortex, phosphoDVP main experiment (Fig. 5e-f).',
          neuron_src, sites, samples)


if __name__ == '__main__':
    main()
