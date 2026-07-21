# PIPELINE.md — BIOFILM_INSILICO_SUPPORT

## Purpose

This document defines the executable computational pipeline for the current in silico support workflow accompanying the thiosugar–levoglucosenone biofilm manuscript.

The pipeline covers:

- validated ligand inputs
- receptor preparation
- ORCA / xTB quantum-chemical calculations
- CHELPG charge extraction
- ADMET prediction curation
- AutoDock Vina docking
- orthosteric-pose validation
- PLIP interaction profiling
- manuscript-ready tables and figures
- logs, software versions, and checksums

This pipeline explicitly excludes:

- GROMACS molecular dynamics
- MM-GBSA
- trajectory analysis
- wetlab raw-data processing
- microscopy raw-image archives

Molecular dynamics is reserved for a separate dedicated repository and manuscript.

---

## Source-of-truth workspace

Primary repository on the workstation:

```text
/data/BIOFILM_INSILICO_SUPPORT/
```

Canonical input folders:

```text
/data/BIOFILM_INSILICO_SUPPORT/00_validated_inputs/ligands/
/data/BIOFILM_INSILICO_SUPPORT/00_validated_inputs/receptors/
```

Do not use legacy MD directories as public workflow inputs. Legacy/local directories may be used only as clearly marked local references.

---

## Compound identifiers

Use only the following compound identifiers:

| Compound | Repository name | Formula | MW |
|---|---|---:|---:|
| 8 | compound_08 | C33H45N3O19S2 | 851.85 |
| 9 | compound_09 | C18H25N3O9S2 | 491.53 |
| 11 | compound_11 | C8H9N3O3S2 | 259.30 |

Do not use obsolete local labels such as 16a, 20, or 23 in public-facing files.

---

## Receptor identifiers

Use only the following receptor identifiers:

| Target | PDB ID | Repository name |
|---|---|---|
| LasR | 2UV0 | LasR_2UV0 |
| PqsR / MvfR | 4JVI | PqsR_4JVI |

---

## Pipeline overview

```text
validated SMILES
      |
      v
SMILES validation and checksum
      |
      v
3D structure generation / geometry verification
      |
      v
ORCA / xTB calculations
      |
      v
CHELPG charge extraction
      |
      v
ADMET raw exports and processed summaries
      |
      v
receptor and ligand PDBQT preparation
      |
      v
AutoDock Vina docking
      |
      v
orthosteric validation by reference-ligand distance
      |
      v
PLIP interaction profiling
      |
      v
manuscript tables, figures, methods text
      |
      v
software-version freeze and checksums
```

---

## Step 0 — Repository sanity check

### Input

Existing directory:

```text
/data/BIOFILM_INSILICO_SUPPORT/
```

Expected top-level structure:

```text
00_validated_inputs/
01_qm_orca_xtb/
02_admet/
03_docking/
04_scripts/
05_manuscript_outputs/
99_logs_versions_checksums/
```

### Output

No new data output. This step verifies that the workspace is structurally complete.

### Quality checkpoint

The repository must not contain GROMACS production files, trajectories, MM-GBSA outputs, or obsolete local compound labels.

---

## Step 1 — Validate ligand SMILES

### Input

```text
00_validated_inputs/ligands/compound_08/compound_08.smiles
00_validated_inputs/ligands/compound_09/compound_09.smiles
00_validated_inputs/ligands/compound_11/compound_11.smiles
```

### Expected current SMILES

Compound 08:

```text
S=C(N)N/N=C1C[C@H](S[C@@H]2[C@H](OC(C)=O)[C@@H](OC(C)=O)[C@H](O[C@@H]3[C@H](OC(C)=O)[C@@H](OC(C)=O)[C@H](OC(C)=O)[C@@H](COC(C)=O)O3)[C@@H](COC(C)=O)O2)[C@@H]4O[C@H]/1OC4
```

Compound 09:

```text
S=C(N)N/N=C1C[C@H](S[C@@H]2[C@H](OC(C)=O)[C@@H](OC(C)=O)[C@H](OC(C)=O)CO2)[C@@H]3O[C@H]/1OC3
```

Compound 11:

```text
O=C1C[C@H](N(N=C(N)S2)C2=S)[C@@H]3O[C@H]1OC3
```

### Output

```text
99_logs_versions_checksums/checksums/sha256_ligand_smiles.txt
05_manuscript_outputs/tables/table_validated_ligands.csv
```

### Quality checkpoint

Validate that the parsed formula and molecular weight match:

| Compound | Required formula | Required MW |
|---|---|---:|
| 8 | C33H45N3O19S2 | 851.85 |
| 9 | C18H25N3O9S2 | 491.53 |
| 11 | C8H9N3O3S2 | 259.30 |

For compounds 8 and 9, confirm that the thiosemicarbazone motif is present:

```text
C=N–NH–C(=S)–NH2
```

Do not use N2 formulas for compounds 8 or 9.

---

## Step 2 — Generate and verify 3D ligand geometries

### Input

Validated SMILES from Step 1.

### Output

```text
01_qm_orca_xtb/compound_08/geometries/
01_qm_orca_xtb/compound_09/geometries/
01_qm_orca_xtb/compound_11/geometries/
```

Recommended file types:

```text
*.sdf
*.mol
*.xyz
```

### Quality checkpoint

Inspect the generated 3D structures before quantum-chemical calculations.

Check:

- stereochemistry of carbohydrate centers
- levoglucosenone stereocentres
- hydrazone C=N geometry
- thiosemicarbazone N–N connectivity in compounds 8 and 9
- absence of broken rings or chemically implausible valence states
- no accidental conversion of compound 9 into a galactose derivative

SMILES encodes stereochemical configuration, not a frozen carbohydrate ring conformation. The final 3D conformation must be checked after embedding and after optimization.

---

## Step 3 — ORCA / xTB geometry and electronic-structure calculations

### Input

3D ligand geometries from Step 2.

### Compound-specific workflow

| Compound | Geometry / electronic workflow | Interpretation rule |
|---|---|---|
| compound_08 | GFN2-xTB pre-optimization followed by omegaB97X-D3/def2-SVP/CPCM(water) optimization; Cartesian-Hessian restart; analytical frequency calculation | converged local minimum with 0 imaginary frequencies; directly comparable with 09/11 at the finalized DFT level |
| compound_09 | omegaB97X-D3/def2-SVP/CPCM(water) optimization and analytical frequency calculation | converged local minimum with 0 imaginary frequencies; directly comparable with 08/11 |
| compound_11 | omegaB97X-D3/def2-SVP/CPCM(water) optimization and analytical frequency calculation | converged local minimum with 0 imaginary frequencies; directly comparable with 08/09 |

The current validated HOMO-LUMO gaps are 9.3309 eV for compound_08, 9.4330 eV for compound_09, and 8.9703 eV for compound_11. This ordering is not a ranking of biological activity, toxicity, or antibiofilm potency.

This section was updated after completion and QC of STEPS 3–4 and supersedes the pre-computation workflow description.

### Input folders

```text
01_qm_orca_xtb/compound_08/inputs/
01_qm_orca_xtb/compound_09/inputs/
01_qm_orca_xtb/compound_11/inputs/
```

### Output folders

```text
01_qm_orca_xtb/compound_08/outputs/
01_qm_orca_xtb/compound_09/outputs/
01_qm_orca_xtb/compound_11/outputs/
```

### Quality checkpoint

Record for every calculation:

- input file
- charge and multiplicity
- method
- basis set
- solvent model
- ORCA or xTB version
- number of CPU cores
- memory settings
- start date and completion date
- convergence status
- imaginary-frequency status, if relevant
- output checksum

---

## Step 4 — CHELPG charge calculation and extraction

### Input

Optimized or finalized geometries from Step 3.

### Charge workflow

Use a consistent charge level for all compounds:

```text
HF/6-31G*
```

### Output

```text
01_qm_orca_xtb/compound_08/charges/
01_qm_orca_xtb/compound_09/charges/
01_qm_orca_xtb/compound_11/charges/
05_manuscript_outputs/tables/table_chelpg_charges.csv
05_manuscript_outputs/tables/table_pharmacophoric_heteroatom_charges.csv
```

### Quality checkpoint

For compound 11, validate atom identity against the RDKit/ORCA mapping before interpreting CHELPG charges. For the validated structure and current HF/6-31G* gas-phase CHELPG calculation:

```text
N-4 (ORCA atom 4) = +0.31321861 e
N-3 (ORCA atom 5) = -0.45537399 e
exocyclic amino N (ORCA atom 7) = -0.95315983 e
```

The historical value `N-4 ≈ -1.087 e` was not reproduced. Its provenance is unresolved; do not use it as a QC target or current result, and do not retune the method to reproduce it.

Do not revive the obsolete Tyr258-centered hypothesis for compound 11. The current PqsR interaction model for compound 11 involves Leu197 / Ile236 backbone and Gln194.

---

## Step 5 — ADMET raw export collection

### Input

Validated SMILES from Step 1.

### Platforms

```text
SwissADME
ADMETlab 3.0
ProTox-3.0
```

pkCSM is not part of the current correction workflow unless explicitly requested for historical comparison.

### Raw-output folders

```text
02_admet/swissadme/raw/
02_admet/admetlab3/raw/
02_admet/protox/raw/
```

### Required files

```text
swissadme_08.csv
swissadme_09.csv
swissadme_11.csv

ADMETlab3_result_08.csv
ADMETlab3_result_09.csv
ADMETlab3_result_11.csv

ProTox-3.0 - Prediction of TOXicity of chemicals_08.csv
ProTox-3.0 - Prediction of TOXicity of chemicals_09.csv
ProTox-3.0 - Prediction of TOXicity of chemicals_11.csv
```

### Quality checkpoint

For every webserver result, record:

- platform name
- access date
- input SMILES
- raw export filename
- whether the export is complete
- checksum

Do not use incomplete ProTox PDF exports with unresolved `Calculating...` fields as primary data. Prefer complete CSV exports.

---

## Step 6 — ADMET processing and summary tables

### Input

Raw ADMET exports from Step 5.

### Output

```text
02_admet/swissadme/processed/
02_admet/admetlab3/processed/
02_admet/protox/processed/
02_admet/summary_tables/
05_manuscript_outputs/tables/table_admet_summary.csv
```

### Required summary fields

At minimum, the manuscript ADMET summary should include:

- formula
- molecular weight
- TPSA
- HBA
- HBD
- rotatable bonds
- consensus logP or clearly specified logP model
- solubility
- GI absorption
- BBB permeability
- P-gp status
- CYP liabilities
- ProTox LD50 and toxicity class
- key ProTox endpoint flags
- ADMETlab3 safety and pharmacokinetic flags

### Quality checkpoint

Use consistent source attribution per field. Do not mix SwissADME, ADMETlab3, ProTox, and historical pkCSM values without explicit labeling.

---

## Step 7 — Receptor preparation

### Input receptors

```text
00_validated_inputs/receptors/LasR_2UV0/
00_validated_inputs/receptors/PqsR_4JVI/
```

### Targets

| Target | PDB ID | Reference ligand |
|---|---|---|
| LasR | 2UV0 | OHN or co-crystallized LasR ligand |
| PqsR / MvfR | 4JVI | QZN |

### Output

```text
03_docking/receptors_pdbqt/LasR_2UV0.pdbqt
03_docking/receptors_pdbqt/PqsR_4JVI.pdbqt
```

### Quality checkpoint

Record receptor-preparation decisions:

- waters removed or retained
- co-crystallized ligand removed for docking
- hydrogens added
- charge model
- chain selection
- binding-site center
- reference ligand coordinates retained separately for orthosteric validation

---

## Step 8 — Ligand PDBQT preparation

### Input

Final 3D geometries from quantum-chemical workflow.

### Output

```text
03_docking/ligands_pdbqt/compound_08.pdbqt
03_docking/ligands_pdbqt/compound_09.pdbqt
03_docking/ligands_pdbqt/compound_11.pdbqt
```

### Quality checkpoint

Record:

- source geometry
- conversion tool and version
- protonation assumptions
- Gasteiger charge assignment
- number of rotatable bonds / TORSDOF
- any locked bonds or manually fixed torsions

The PDBQT files must be traceable to the validated ligand structure and documented 3D geometry.

---

## Step 9 — AutoDock Vina docking

### Input

Prepared receptor and ligand PDBQT files.

### Output

```text
03_docking/vina_configs/
03_docking/vina_outputs/
03_docking/summary_tables/table_vina_raw_scores.csv
```

### Recommended docking matrix

| Compound | LasR_2UV0 | PqsR_4JVI |
|---|---|---|
| compound_08 | run | run |
| compound_09 | run | run |
| compound_11 | run | run |

### Quality checkpoint

For every docking run, record:

- receptor
- ligand
- box center
- box dimensions
- exhaustiveness
- seed, if set
- top Vina score
- output pose file
- log file
- command used

Affinity alone is not a valid mechanistic conclusion.

---

## Step 10 — Orthosteric-pose validation

### Input

Docking poses from Step 9 and reference ligand positions from receptor structures.

### Output

```text
03_docking/orthosteric_validation/
05_manuscript_outputs/tables/table_orthosteric_validation.csv
```

### Rule

A docking pose is classified as orthosteric only if:

```text
ligand centre-of-mass distance to reference ligand <= 5 Å
```

If the distance is greater than 5 Å, classify the pose as:

```text
non-specific / surface binding
```

even if the Vina score is favorable.

### Required output fields

```text
compound
receptor
vina_score_kcal_mol
reference_ligand
com_distance_angstrom
classification
pose_file
comment
```

---

## Step 11 — PLIP interaction profiling

### Input

Validated orthosteric poses and selected non-specific poses, if biologically relevant.

### Output

```text
03_docking/plip_outputs/
05_manuscript_outputs/tables/table_plip_interactions.csv
```

### Required fields

```text
compound
receptor
pose_id
classification
residue
chain
interaction_type
distance_angstrom
angle_degrees_if_available
comment
```

### Quality checkpoint

Interpret PLIP data only in the context of orthosteric validation. A PLIP contact from a surface pose should not be used as evidence of receptor antagonism.

---

## Step 12 — Integrated docking interpretation

### Expected current interpretation logic

| Compound | Expected interpretation |
|---|---|
| compound_08 | Matrix-directed activity with validated PqsR pose but no strong evidence of global QS silencing |
| compound_09 | Matrix-directed / defective-matrix phenotype; high raw Vina scores do not count if poses are non-orthosteric |
| compound_11 | Validated orthosteric LasR and PqsR binding; QS repression hypothesis supported computationally but still not causal proof |

### Quality checkpoint

Separate facts from interpretation.

Facts:

- Vina score
- distance to reference ligand
- PLIP contacts
- ADMET prediction
- wetlab phenotype

Interpretations:

- geometric pre-organization
- defective matrix
- pqsA futile cycle
- receptor-independent matrix intercalation

---

## Step 13 — Manuscript tables and figures

### Output folders

```text
05_manuscript_outputs/tables/
05_manuscript_outputs/figures/
05_manuscript_outputs/methods_text/
```

### Expected tables

```text
table_validated_ligands.csv
table_admet_summary.csv
table_qm_fmo_summary.csv
table_chelpg_charges.csv
table_pharmacophoric_heteroatom_charges.csv
table_vina_raw_scores.csv
table_orthosteric_validation.csv
table_plip_interactions.csv
table_docking_summary.csv
```

### Expected methods text

```text
qm_methods.md
admet_methods.md
docking_methods.md
reproducibility_statement.md
```

### Quality checkpoint

Every numerical value in manuscript-facing tables must have a traceable source file.

---

## Step 14 — Software versions and environment freeze

### Output

```text
99_logs_versions_checksums/software_versions/
```

### Required version records

```text
orca_version.txt
xtb_version.txt
openbabel_version.txt
vina_version.txt
plip_version.txt
python_environment.yml
r_environment.txt
system_info.txt
```

### Quality checkpoint

Do not guess software versions. Record exact local outputs.

---

## Step 15 — Checksums and release snapshot

### Output

```text
99_logs_versions_checksums/checksums/
```

### Required checksum files

```text
sha256_ligand_smiles.txt
sha256_receptors.txt
sha256_scripts.txt
sha256_outputs.txt
```

### Recommended command pattern

```bash
sha256sum 00_validated_inputs/ligands/*/* > 99_logs_versions_checksums/checksums/sha256_ligand_smiles.txt
sha256sum 00_validated_inputs/receptors/*/* > 99_logs_versions_checksums/checksums/sha256_receptors.txt
find 04_scripts -type f -exec sha256sum {} \; > 99_logs_versions_checksums/checksums/sha256_scripts.txt
find 05_manuscript_outputs -type f -exec sha256sum {} \; > 99_logs_versions_checksums/checksums/sha256_outputs.txt
```

### Quality checkpoint

Run checksums:

- after validating ligand inputs
- after final receptor preparation
- after final scripts
- before manuscript submission
- before public repository release

---

## Terminology guardrails

Use:

- levoglucosenone
- 6,8-dioxabicyclo[3.2.1]octane
- thiosemicarbazone for compounds 8 and 9
- thioether bridge for the sugar linkage in compounds 8 and 9
- xylopyranose for compound 9
- 1,3,4-thiadiazole aza-adduct for compound 11
- LasR 2UV0
- PqsR / MvfR 4JVI

Do not use:

- thiourea for compounds 8 and 9
- chromene or coumarin for the levoglucosenone scaffold
- galactose for compound 9
- 3IX3 for LasR
- 4JVC for PqsR
- obsolete labels 16a, 20, or 23 in public-facing files

---

## Final release checklist

Before using results in the manuscript or repository release, confirm:

- [ ] final SMILES are present and checksummed
- [ ] formula and MW are verified for compounds 8, 9, and 11
- [ ] ORCA / xTB logs are archived
- [ ] CHELPG charge tables are traceable
- [ ] ADMET raw exports are archived
- [ ] ADMET processed summaries are traceable
- [ ] receptor preparation is documented
- [ ] ligand PDBQT preparation is documented
- [ ] Vina configuration files are archived
- [ ] orthosteric validation table is complete
- [ ] PLIP interaction table is complete
- [ ] manuscript tables are regenerated from current inputs
- [ ] software versions are frozen
- [ ] checksums are updated
- [ ] no GROMACS / MD / MM-GBSA files are included in this repository
