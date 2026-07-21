# STEP 8 ligand PDBQT preparation

Date: 2026-07-16

## Method

- Final frozen STEP 3 QM coordinates were transferred atom-for-atom onto audited MOL topology.
- No bond perception from distances and no geometry optimization were performed.
- The validated neutral protomer/tautomer was retained for each compound.
- Meeko assigned Gasteiger charges and merged standard nonpolar hydrogens.
- REMARK INDEX MAP and deterministic atom renaming were retained.
- Chemistry-aware SMARTS rules rigidified conjugated exocyclic bonds misclassified by default Meeko.

## Results

| Compound | Formula | Input atoms | PDBQT atoms | TORSDOF | Gasteiger charge sum | SDF max coordinate difference (Å) | PDBQT max coordinate difference (Å) |
|---|---:|---:|---:|---:|---:|---:|---:|
| compound_08 | C33H45N3O19S2 | 102 | 60 | 20 | 0.003000 | 0.000000 | 0.000495 |
| compound_09 | C18H25N3O9S2 | 57 | 35 | 8 | 0.005000 | 0.000000 | 0.000500 |
| compound_11 | C8H9N3O3S2 | 25 | 18 | 1 | 0.000000 | 0.000000 | 0.000495 |

## Torsion rules

### compound_08

- SMARTS: `[C;X3](=[N;X2]-[N;X3]-[C;X3](=[S;X1])-[N;X3])`
- Bond indices: `2 3`
- Rule: Rigidify conjugated hydrazone N-N bond in C=N-NH-C(=S)-NH2.
- Final TORSDOF: 20

### compound_09

- SMARTS: `[C;X3](=[N;X2]-[N;X3]-[C;X3](=[S;X1])-[N;X3])`
- Bond indices: `2 3`
- Rule: Rigidify conjugated hydrazone N-N bond in C=N-NH-C(=S)-NH2.
- Final TORSDOF: 8

### compound_11

- SMARTS: `[N;X3;H2]-[c;R]`
- Bond indices: `1 2`
- Rule: Rigidify conjugated exocyclic NH2-thiadiazole C-N bond; retain ring-N-scaffold sigma rotor.
- Final TORSDOF: 1

## QC status

All compounds passed formula, atom-order, stereochemical-identity, charge, TORSDOF, atom-count, and coordinate-preservation checks.

Production docking was not performed in this step.

## Reproducibility scripts

| Script | Role | SHA-256 |
|---|---|---|
| `04_scripts/python/step8b_audit_ligand_topology_mapping.py` | Topology and atom-order audit | `6e5ddf7c344c1b8ff3f39e7d91ff85b2749c0521522861c97996814607694bf2` |
| `04_scripts/python/step8b_map_meeko_rotors.py` | Default Meeko rotor mapping | `713b8c3088ba7b0be790d20bb98b676f2bc45d94fdeb1a28194135e5adc46535` |
| `04_scripts/python/step8b_validate_proposed_torsion_rules.py` | Validation of chemistry-aware torsion rules | `a34be9b636fc75d206d46c24d5f83ccf371eb1fd3205e2a159278c667b5d0107` |
| `04_scripts/python/step8_prepare_final_qm_ligands.py` | Canonical ligand SDF/PDBQT preparation | `8203a008e3a86fd7a38219094b1b778c4137b254f304d4e6cc93e49a17b7501f` |
| `04_scripts/python/step8_finalize_ligand_documentation.py` | Documentation and checksum finalization | `7d6cc08c51e37b3e8c9c2b6aae588117feef448e54cccb5523e33457dc609650` |

The preparation and audit scripts were archived in the repository before STEP 8 backup.
The STEP 8 checksum manifest was regenerated after adding these scripts and updating this documentation.
