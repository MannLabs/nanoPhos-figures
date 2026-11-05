"""Helper to dump each figure panel's exact plotted source data into the PRIDE
MetaInfo workbook (one sheet per panel), matching MetaInfo_figures_checklist_v01.xlsx.

Usage in a figure notebook (run in Jupyter, after the panel's plotting dataframe is built):

    from metainfo_export import dump_panel
    dump_panel(df, "Figure 2b")          # df = the exact dataframe the plot consumes
    dump_panel(df, "Figure 2f", notes={"x-axis": "PC1 (0.73)"})   # values not in df

The first call creates MetaInfo_figures_checklist_v05.xlsx; subsequent calls append /
replace that panel's sheet, so notebooks can be run in any order. Sheet names must
match the figure panel labels.
"""
from __future__ import annotations
import os
import pandas as pd

DEFAULT_PATH = os.path.join("source_data", "MetaInfo_figures_checklist_v05.xlsx")


def dump_panel(df: pd.DataFrame, sheet: str, path: str = DEFAULT_PATH,
               index: bool = False, notes: dict | None = None) -> None:
    """Write `df` to `sheet` in the MetaInfo v05 workbook (create or replace the sheet).

    `notes` holds plotted values that are not columns of `df` (axis titles with explained
    variance, test statistics, fit parameters, thresholds, ordering); they are written as
    a two-column Note / Value block one column to the right of the data.
    """
    df = pd.DataFrame(df).copy()
    # openpyxl sheet-name limit is 31 chars
    sheet = sheet[:31]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    exists = os.path.exists(path)
    with pd.ExcelWriter(path, engine="openpyxl", mode="a" if exists else "w",
                        **({"if_sheet_exists": "replace"} if exists else {})) as xl:
        df.to_excel(xl, sheet_name=sheet, index=index)
        if notes:
            ws, col = xl.sheets[sheet], df.shape[1] + int(index) + 2
            ws.cell(row=1, column=col, value="Note")
            ws.cell(row=1, column=col + 1, value="Value")
            for r, (k, v) in enumerate(notes.items(), start=2):
                ws.cell(row=r, column=col, value=k)
                ws.cell(row=r, column=col + 1, value=v)
    print(f"  [MetaInfo] wrote '{sheet}'  ({df.shape[0]} rows x {df.shape[1]} cols"
          f"{f', {len(notes)} notes' if notes else ''})")
