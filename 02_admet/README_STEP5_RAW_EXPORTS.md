# Step 5 — ADMET raw export provenance

Access date for all exports: **2026-07-02**.

## Tools

- SwissADME — web-service version/build not stated in the CSV export.
- ADMETlab 3.0.
- ProTox 3.0.

## Structural verification

SwissADME and ADMETlab exports contain SMILES fields. For compounds 08, 09 and 11, these structures were canonicalized with RDKit and matched the current validated isomeric SMILES.

ProTox 3.0 CSV exports do not contain the submitted SMILES. Their compound assignment is therefore documented by provenance rather than by direct structural reconstruction from the CSV.

## Molecular-mass convention

SwissADME reports average molecular weight. ADMETlab 3.0 reports values consistent with monoisotopic mass. Raw values are preserved unchanged. Manuscript molecular weights remain:

- compound_08: 851.85 g/mol
- compound_09: 491.53 g/mol
- compound_11: 259.30 g/mol

## Data-integrity rule

Raw exports are immutable source files. Any harmonized tables must be generated separately and must retain source-tool attribution.

## Export completeness

All nine raw exports were checked for completeness.

- SwissADME: one complete compound record per CSV.
- ADMETlab 3.0: one complete compound record per CSV.
- ProTox 3.0: 45 prediction rows per CSV and no unresolved `Calculating...` fields.

The literal submitted SMILES, raw export filename, completeness status, verification basis, and checksum are recorded in `admet_raw_export_provenance.csv`.
