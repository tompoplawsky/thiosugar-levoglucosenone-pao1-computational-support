# Computational environment — public repository version

## Scope

This document describes the computational environment for the current in silico support workflow accompanying the thiosugar–levoglucosenone biofilm manuscript.

Included:

- ORCA / xTB quantum-chemical calculations
- CHELPG charge analysis
- ADMET prediction processing
- AutoDock Vina docking
- PLIP interaction profiling
- scripts, tables, figures, logs, checksums

Excluded:

- GROMACS molecular dynamics
- MM-GBSA
- MDAnalysis trajectory analysis
- wetlab source data
- local legacy MD workspaces

Molecular dynamics is reserved for a separate dedicated repository and manuscript.

## Hardware environment

| Component | Specification |
|---|---|
| Workstation | Dell Precision 5820 Tower |
| CPU | Intel Xeon W-2145, 8 cores / 16 threads |
| RAM | 64 GB DDR4 ECC |
| GPU | NVIDIA GeForce RTX 3060, 12 GB VRAM |
| Operating system | Linux Mint 22.3, based on Ubuntu 24.04 LTS |

The GPU is not a required component of the current ORCA/docking/ADMET support workflow.

## Software records

Store software metadata under:

`99_logs_versions_checksums/software_versions/`

Required version records:

| Tool | Purpose | Version file |
|---|---|---|
| ORCA | DFT, HF single points, CHELPG charges | orca_version.txt |
| xTB | GFN2-xTB pre-optimization | xtb_version.txt |
| Open Babel | structure conversion | openbabel_version.txt |
| AutoDock Vina | docking | vina_version.txt |
| PLIP | interaction profiling | plip_version.txt |
| Python | parsing, summaries, plotting | python_environment.yml |
| R | optional table/figure/statistical utilities | r_environment.txt |
| OS | workstation metadata | system_info.txt |

## ORCA / xTB workflow

| Compound | Geometry workflow | Interpretation note |
|---|---|---|
| compound_08 | GFN2-xTB pre-optimization + HF/6-31G* single point | HOMO-LUMO gap not directly comparable with compounds 09 and 11 |
| compound_09 | omegaB97X-D3/def2-SVP/CPCM(water) | comparable with compound 11 |
| compound_11 | omegaB97X-D3/def2-SVP/CPCM(water) | comparable with compound 09 |

CHELPG charges should be calculated consistently at HF/6-31G* for all three compounds.

## ADMET workflow

Supported platforms:

- SwissADME
- ADMETlab 3.0
- ProTox-3.0

Rules:

- keep raw exports in `raw/`
- keep cleaned files in `processed/`
- record access date and tool version or platform name
- date-stamp re-runs
- do not overwrite raw exports

## Docking workflow

Targets:

| Target | PDB ID |
|---|---|
| LasR | 2UV0 |
| PqsR / MvfR | 4JVI |

A docking pose is classified as orthosteric only if the ligand centre of mass is within 5 angstroms of the co-crystallized reference ligand.

Affinity alone is insufficient. Surface poses should be reported as non-specific even if the raw Vina score is favorable.

## Repository boundary

This repository supports the present wetlab manuscript by providing reproducible computational evidence from ORCA/xTB, ADMET, docking, and PLIP analyses.

It is not the repository for the planned molecular-dynamics manuscript.
