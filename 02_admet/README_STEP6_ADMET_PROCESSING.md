# Step 6 — ADMET processing and integrated summary

Access date for all web-server predictions: **2026-07-02**.

## Source allocation

- Validated project inputs: molecular formula and manuscript molecular weight.
- SwissADME: TPSA, HBA, HBD, rotatable bonds, consensus logP, ESOL solubility, GI absorption, BBB permeability, P-gp substrate status, and qualitative CYP calls.
- ADMETlab 3.0: regression outputs and endpoint-category probabilities for pharmacokinetic and safety models.
- ProTox 3.0 CSV: endpoint predictions and probabilities.
- ProTox 3.0 PDF, page 1: predicted oral LD50, toxicity class, average similarity, and prediction accuracy.

## Interpretation rule

ADMETlab numerical outputs are retained as exported. Classification-model values are recorded as probabilities of the model-specific category 1. No universal binary threshold was applied across all endpoints because category definitions differ between models.

SwissADME values marked `n/d` remain `n/d`; they were not replaced with predictions from another platform.

ProTox acute oral toxicity values were visually verified from rendered page 1 of each PDF because this graphical block is not present in the extractable PDF text layer.

## Limitations

All ADMET and toxicity results are computational predictions and do not replace experimental pharmacokinetic or cytotoxicity measurements. The ProTox acute oral predictions have prediction accuracy of 23% and average similarity of approximately 37–40%, so toxicity classes must be interpreted cautiously.

No pkCSM values were incorporated.

## Manuscript and detailed tables

`02_admet/summary_tables/admet_integrated_summary.csv` contains the complete 65-column integrated dataset.

`05_manuscript_outputs/tables/table_admet_detailed.csv` is an identical detailed copy intended for Supporting Information.

`05_manuscript_outputs/tables/table_admet_summary.csv` is the condensed manuscript table. ADMETlab probabilities are rounded to three decimal places only in this condensed presentation; unrounded values remain available in the detailed tables. No universal binary threshold was applied.

ProTox endpoint calls retain the exported Active/Inactive classification together with the corresponding probability.
