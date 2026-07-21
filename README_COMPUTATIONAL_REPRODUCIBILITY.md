# Computational reproducibility checklist

## Repository boundary

This repository covers only ORCA/xTB, ADMET, docking, PLIP, scripts, and manuscript computational outputs for the current wetlab manuscript.

It excludes GROMACS/MD and MM-GBSA, which are reserved for a separate MD manuscript.

## Minimal reproducibility requirements

For every computational result, record:

- input file
- software version
- command or webserver settings
- date
- raw output
- processed output
- checksum

## Required checksums

Recommended checksum files:

- `sha256_ligand_smiles.txt`
- `sha256_receptors.txt`
- `sha256_scripts.txt`
- `sha256_outputs.txt`

## Raw outputs

Do not overwrite raw outputs. Re-runs must be date-stamped.

## Public exclusions

Do not include:

- trajectory files
- GROMACS production files
- MM-GBSA files
- raw wetlab archives
- obsolete local compound labels
- user-specific desktop paths
