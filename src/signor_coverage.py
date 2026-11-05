"""Supplementary Table 2: coverage of the curated SIGNOR EGFR-signaling pathway.

Pathway sites: unique (protein, residue, position) phosphorylation / dephosphorylation sites
on proteins in the SIGNOR "EGFR Signaling" export (data/SIGNOR-EGF_25_06_26.tsv, downloaded
25 June 2026); protein families and complexes are excluded. A site is covered at an input
when it is a Class I site (per-(site, run) localization >= 0.75) in at least one run of the
EGF dilution series at that input (3 EGF-stimulated + 3 unstimulated runs, per-input searches).

Run from the repo root:  python src/signor_coverage.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from core import process_ptm_site_report

CLASS_I_CUTOFF = 0.75
INPUTS = [10, 20, 50, 100, 200, 500, 1000]
SIGNOR = Path('data/SIGNOR-EGF_25_06_26.tsv')
EGF_REPORTS = Path('pride_data/analysis_data/revision/figure2')
OUT = Path('tables/Supplementary_Table_2.xlsx')
PROC_META = {'Protein_group', 'Gene_group', 'PTM_0_aa', 'PTM_pos', 'PTM_mult123',
             'PTM_flank', 'PTM_Collapse_key', 'PTM_localization', 'UPD_seq'}
AA = {'Ser': 'S', 'Thr': 'T', 'Tyr': 'Y'}


def pathway_sites():
    sg = pd.read_csv(SIGNOR, sep='\t')
    sg = sg[sg['TYPEB'].eq('protein') & sg['MECHANISM'].isin(['phosphorylation', 'dephosphorylation'])
            & sg['RESIDUE'].notna()]
    sg = sg.assign(aa=sg['RESIDUE'].str[:3].map(AA), pos=sg['RESIDUE'].str[3:].astype(int))
    kinases = sg.groupby(['IDB', 'aa', 'pos'])['ENTITYA'].agg(lambda s: ';'.join(sorted(set(s))))
    sites = sg.drop_duplicates(['IDB', 'aa', 'pos'])[['ENTITYB', 'IDB', 'aa', 'pos']]
    sites = sites.assign(kinases=[kinases[k] for k in zip(sites['IDB'], sites['aa'], sites['pos'])])
    return sites.sort_values(['ENTITYB', 'pos']).reset_index(drop=True)


def covered_sites(report, keys):
    """Pathway sites that are Class I in at least one run of `report`."""
    sd = process_ptm_site_report(pd.read_csv(report, sep='\t', low_memory=False),
                                 cutoff=CLASS_I_CUTOFF)['site_data']
    runs = [c for c in sd.columns if c not in PROC_META]
    sd = sd[sd[runs].notna().any(axis=1)]
    return {(acc, aa, int(pos)) for pg, aa, pos in zip(sd['Protein_group'], sd['PTM_0_aa'], sd['PTM_pos'])
            for acc in str(pg).split(';') if (acc, aa, int(pos)) in keys}


def main():
    sites = pathway_sites()
    keys = set(zip(sites['IDB'], sites['aa'], sites['pos']))
    covered = {ng: covered_sites(next(EGF_REPORTS.glob(f'*withEGF*_{ng}ng_Report.tsv')), keys) for ng in INPUTS}

    matrix = pd.DataFrame({'substrate': sites['ENTITYB'], 'acc': sites['IDB'],
                           'site': sites['aa'] + sites['pos'].astype(str), 'aa': sites['aa'],
                           'kinases': sites['kinases']})
    for ng in INPUTS:
        matrix[f'{ng}ng'] = ['yes' if k in covered[ng] else 'no'
                             for k in zip(sites['IDB'], sites['aa'], sites['pos'])]
    detected = matrix[[f'{ng}ng' for ng in INPUTS]].eq('yes')
    matrix['n_inputs_detected'] = detected.sum(axis=1)
    matrix['min_input_detected_ng'] = detected.apply(lambda r: next((ng for ng, d in zip(INPUTS, r) if d), None), axis=1)

    union = set().union(*covered.values())
    summary = pd.DataFrame([{'input_ng': ng, 'n_sites': len(keys), 'nanoPhos_covered': len(covered[ng]),
                             'nanoPhos_pct': round(100 * len(covered[ng]) / len(keys), 1)} for ng in INPUTS]
                           + [{'input_ng': 'union(all)', 'n_sites': len(keys), 'nanoPhos_covered': len(union),
                               'nanoPhos_pct': round(100 * len(union) / len(keys), 1)}])
    with pd.ExcelWriter(OUT) as xl:
        summary.to_excel(xl, sheet_name='Summary', index=False)
        matrix.to_excel(xl, sheet_name='Sites_across_dilutions', index=False)
    print(summary.to_string(index=False))
    print(f'{OUT}: {len(keys)} pathway sites, {len(union)} covered across the dilution series')


if __name__ == '__main__':
    main()
