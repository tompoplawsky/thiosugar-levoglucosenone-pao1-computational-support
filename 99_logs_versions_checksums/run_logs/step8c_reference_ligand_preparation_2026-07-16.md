# STEP 8C reference-ligand preparation

Date: 2026-07-16

## Method

- Deposited heavy-atom coordinates were preserved exactly.
- Bond topology was reconstructed from embedded RCSB chemical-component dictionaries.
- No distance-based bond perception or geometry optimization was used.
- Hydrogens were added by RDKit without moving heavy atoms.
- Neutral deposited chemical-component states were retained without enumeration.
- Meeko assigned Gasteiger charges and merged standard nonpolar hydrogens.
- OHN used default Meeko torsions.
- QZN used a validated SMARTS rule that rigidified only the conjugated exocyclic amino N–N bond.
- QZN site B corresponds to biological-assembly chain A-2 mapped to working PDB chain B.

## Results

| Reference | Formula | Heavy atoms | Added H | PDBQT atoms | TORSDOF | Rounded charge sum | Source Δmax (Å) | SDF Δmax (Å) | PDBQT Δmax (Å) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| LasR_2UV0_OHN_chainF | C16H27NO4 | 21 | 27 | 22 | 11 | -0.000000 | 0.000000 | 0.000000 | 0.000000 |
| PqsR_4JVI_QZN_siteA | C17H24ClN3O | 22 | 24 | 24 | 8 | -0.001000 | 0.000000 | 0.000000 | 0.000000 |
| PqsR_4JVI_QZN_siteB | C17H24ClN3O | 22 | 24 | 24 | 8 | -0.001000 | 0.000000 | 0.000000 | 0.000000 |

## Reference-specific details

### LasR_2UV0_OHN_chainF

- Label: LasR 2UV0 OHN chain F
- Compound ID: OHN
- Coordinate source: `00_validated_inputs/receptors/LasR_2UV0/LasR_2UV0_OHN_chainF_reference.pdb`
- Coordinate mmCIF chain: `F`
- Working PDB chain: `F`
- Isomeric SMILES: `CCCCCCCCCC(=O)CC(=O)N[C@H]1CCOC1=O`
- Chiral centers: `[(15, 'S')]`
- Torsion rule: Default Meeko torsions; the conjugated amide C(=O)-N bond is automatically rigid.
- Final TORSDOF: 11

### PqsR_4JVI_QZN_siteA

- Label: PqsR 4JVI QZN site A
- Compound ID: QZN
- Coordinate source: `00_validated_inputs/receptors/PqsR_4JVI/PqsR_4JVI_QZN_chainA_reference.pdb`
- Coordinate mmCIF chain: `A`
- Working PDB chain: `A`
- Isomeric SMILES: `CCCCCCCCCc1nc2cc(Cl)ccc2c(=O)n1N`
- Chiral centers: `[]`
- Torsion rule: Rigidify only the conjugated exocyclic amino N-N bond; retain the alkyl-chain attachment rotor.
- Final TORSDOF: 8

### PqsR_4JVI_QZN_siteB

- Label: PqsR 4JVI QZN site B
- Compound ID: QZN
- Coordinate source: `00_validated_inputs/receptors/PqsR_4JVI/PqsR_4JVI_QZN_chainB_reference.pdb`
- Coordinate mmCIF chain: `A-2`
- Working PDB chain: `B`
- Isomeric SMILES: `CCCCCCCCCc1nc2cc(Cl)ccc2c(=O)n1N`
- Chiral centers: `[]`
- Torsion rule: Rigidify only the conjugated exocyclic amino N-N bond; retain the alkyl-chain attachment rotor.
- Final TORSDOF: 8

## QC status

All reference ligands passed formula, formal-charge, atom-count, topology, coordinate-preservation, PDBQT, charge, and torsion checks.

Redocking and production docking were not performed in STEP 8C.
