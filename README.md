# BIOFILM_INSILICO_SUPPORT

Computational support repository for the current thiosugar–levoglucosenone biofilm manuscript.

This repository contains only the in silico support workflow used to rationalize experimental antibiofilm and quorum-sensing phenotypes of compounds 08, 09, and 11 against *Pseudomonas aeruginosa* PAO1.

## Scope

Included:

- validated ligand structures
- validated receptor structures
- ORCA / xTB quantum-chemical calculations
- CHELPG charge analysis
- ADMET prediction outputs
- AutoDock Vina docking
- PLIP protein-ligand interaction profiling
- scripts
- manuscript-ready computational tables and figures
- software versions, logs, checksums

Excluded:

- wetlab source data
- GROMACS molecular dynamics
- MM-GBSA
- MDAnalysis trajectory analysis
- legacy local MD directories

## Compounds

| Compound | Structural class | Formula | MW |
|---|---|---:|---:|
| compound_08 | peracetylated disaccharide thiosugar-levoglucosenone thiosemicarbazone | C33H45N3O19S2 | 851.85 |
| compound_09 | peracetylated xylopyranose thiosugar-levoglucosenone thiosemicarbazone | C18H25N3O9S2 | 491.53 |
| compound_11 | rigid 1,3,4-thiadiazole aza-adduct | C8H9N3O3S2 | 259.30 |

## Receptors

| Target | PDB ID | Use |
|---|---|---|
| LasR | 2UV0 | orthosteric docking |
| PqsR / MvfR | 4JVI | orthosteric docking |

## Directory overview

BIOFILM_INSILICO_SUPPORT/
  00_validated_inputs/
  01_qm_orca_xtb/
  02_admet/
  03_docking/
  04_scripts/
  05_manuscript_outputs/
  99_logs_versions_checksums/

Detailed directory rules are provided in `REPO_STRUCTURE.md`.

## Docking validation rule

Docking affinity alone is not sufficient.

A pose is considered orthosteric only when the ligand centre of mass is within 5 angstroms of the co-crystallized reference ligand.

Surface poses or poses outside this threshold should be reported as non-specific regardless of raw Vina score.

## Reproducibility policy

All computational outputs should be traceable to:

- validated input structure
- software version
- command or web-tool settings
- date of calculation or prediction
- raw output file
- processed summary table
- checksum

Raw exports should not be overwritten. Re-runs should be date-stamped.
