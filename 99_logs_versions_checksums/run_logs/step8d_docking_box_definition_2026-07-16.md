# STEP 8D docking-box definition

Date: 2026-07-16

## Final box policy

- Box shape: cube.
- Box dimensions: 30.000 × 30.000 × 30.000 Å.
- Box volume: 27,000 Å³.
- Center: crystallographic reference-ligand heavy-atom geometric centroid.
- The same receptor-specific box will be used for reference redocking and later docking of compounds 8, 9, and 11.
- All receptor heavy atoms within 6.0 Å of the crystallographic reference ligand were required to lie inside the box.
- All ligands were required to fit at the box center under arbitrary rotation, using maximum centroid radius.
- Vina search parameters were not defined and docking was not run.

## Validated systems

| System | Center x | Center y | Center z | Size x | Size y | Size z | Pocket atoms | Pocket clearance | Largest ligand | Radius | Rotation clearance |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|
| LasR_2UV0 | 55.258 | 27.138 | 29.990 | 30.000 | 30.000 | 30.000 | 158 | 3.664 Å | compound_08 | 11.022 Å | 3.978 Å |
| PqsR_4JVI | -33.489 | 56.865 | 9.100 | 30.000 | 30.000 | 30.000 | 110 | 2.536 Å | compound_08 | 11.022 Å | 3.978 Å |

## QC status

Both receptor-specific boxes passed center, reference-ligand containment, pocket containment, arbitrary-rotation ligand-fit, and volume checks.

No redocking or production docking was performed.
