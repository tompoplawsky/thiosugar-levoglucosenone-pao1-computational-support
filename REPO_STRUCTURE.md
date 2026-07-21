# Repository structure — BIOFILM_INSILICO_SUPPORT

## Scope

This repository contains the computational support workflow for the current thiosugar–levoglucosenone biofilm manuscript.

Included:

- validated ligand and receptor inputs
- ORCA / xTB quantum-chemical calculations
- CHELPG charge analysis
- ADMET predictions
- AutoDock Vina docking
- PLIP interaction profiling
- scripts
- manuscript-ready computational tables and figures
- logs, versions, checksums

Excluded:

- GROMACS molecular dynamics
- MM-GBSA
- MDAnalysis trajectory analysis
- wetlab source data
- raw microscopy archives
- legacy local MD directories

## Directory structure

BIOFILM_INSILICO_SUPPORT/
  00_validated_inputs/
  01_qm_orca_xtb/
  02_admet/
  03_docking/
  04_scripts/
  05_manuscript_outputs/
  99_logs_versions_checksums/

## 00_validated_inputs

Validated starting structures.

Ligands:

- compound_08
- compound_09
- compound_11

Receptors:

- LasR_2UV0
- PqsR_4JVI

Compound identities:

| Compound | Description | Formula | MW |
|---|---|---:|---:|
| compound_08 | peracetylated disaccharide thiosugar-levoglucosenone thiosemicarbazone | C33H45N3O19S2 | 851.85 |
| compound_09 | peracetylated xylopyranose thiosugar-levoglucosenone thiosemicarbazone | C18H25N3O9S2 | 491.53 |
| compound_11 | rigid 1,3,4-thiadiazole aza-adduct | C8H9N3O3S2 | 259.30 |

## 01_qm_orca_xtb

For each compound:

- inputs
- outputs
- geometries
- charges

## 02_admet

Platforms:

- swissadme
- admetlab3
- protox
- summary_tables

Keep raw and processed files separate.

## 03_docking

Subfolders:

- receptors_pdbqt
- ligands_pdbqt
- vina_configs
- vina_outputs
- plip_outputs
- orthosteric_validation
- summary_tables

Orthosteric classification requires centre-of-mass distance <= 5 angstroms from the co-crystallized reference ligand.

## 04_scripts

Subfolders:

- bash
- python
- r

Scripts should be deterministic and should not contain user-specific absolute paths unless explicitly marked as local helpers.

## 05_manuscript_outputs

Subfolders:

- tables
- figures
- methods_text

## 99_logs_versions_checksums

Subfolders:

- software_versions
- run_logs
- checksums

## Naming policy

Allowed compound names:

- compound_08
- compound_09
- compound_11

Allowed receptor names:

- LasR_2UV0
- PqsR_4JVI

Do not use obsolete local compound labels in public-facing files.
