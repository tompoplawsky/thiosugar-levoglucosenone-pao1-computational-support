#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import importlib.metadata
import platform
import sys
from pathlib import Path

DATE = "2026-07-15"


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_new(path: Path, content: str) -> None:
    if path.exists():
        fail(f"refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def pkg_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not detected"


def count_records(path: Path, prefixes: tuple[str, ...]) -> int:
    count = 0
    with path.open("r", encoding="ascii") as handle:
        for line in handle:
            if line.startswith(prefixes):
                count += 1
    return count


def pdbqt_charge(path: Path) -> float:
    total = 0.0
    with path.open("r", encoding="ascii") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                total += float(line.split()[-2])
    return total


def main() -> None:
    root = Path("/data/BIOFILM_INSILICO_SUPPORT").resolve()
    if not root.is_dir():
        fail(f"project root not found: {root}")

    versions_dir = root / "99_logs_versions_checksums/software_versions"
    checksums_dir = root / "99_logs_versions_checksums/checksums"
    run_logs_dir = root / "99_logs_versions_checksums/run_logs"

    versions_path = versions_dir / f"step7_software_versions_{DATE}.txt"
    report_path = run_logs_dir / f"step7_receptor_preparation_{DATE}.md"
    manifest_path = checksums_dir / f"sha256_step7_receptors_{DATE}.txt"

    for path in (versions_path, report_path, manifest_path):
        if path.exists():
            fail(f"target already exists: {path}")

    relpaths = [
        "00_validated_inputs/receptors/LasR_2UV0/2UV0_rcsb_original.cif",
        "00_validated_inputs/receptors/LasR_2UV0/2UV0_rcsb_original.pdb",
        "00_validated_inputs/receptors/LasR_2UV0/LasR_2UV0_prepared_chainF.pdb",
        "00_validated_inputs/receptors/LasR_2UV0/LasR_2UV0_OHN_chainF_reference.pdb",
        "00_validated_inputs/receptors/LasR_2UV0/LasR_2UV0_HOH2134_chainF_reference.pdb",
        "00_validated_inputs/receptors/PqsR_4JVI/4JVI_rcsb_original.cif",
        "00_validated_inputs/receptors/PqsR_4JVI/4JVI_rcsb_original.pdb",
        "00_validated_inputs/receptors/PqsR_4JVI/4JVI_rcsb_biological_assembly1.cif",
        "00_validated_inputs/receptors/PqsR_4JVI/4JVI_rcsb_biological_assembly1.pdb",
        "00_validated_inputs/receptors/PqsR_4JVI/PqsR_4JVI_prepared_dimer_AB.pdb",
        "00_validated_inputs/receptors/PqsR_4JVI/PqsR_4JVI_QZN_chainA_reference.pdb",
        "00_validated_inputs/receptors/PqsR_4JVI/PqsR_4JVI_QZN_chainB_reference.pdb",
        "00_validated_inputs/receptors/PqsR_4JVI/PqsR_4JVI_HOH501_chainA_reference.pdb",
        "00_validated_inputs/receptors/PqsR_4JVI/PqsR_4JVI_HOH501_chainB_reference.pdb",
        "03_docking/receptors_pdbqt/LasR_2UV0.pdbqt",
        "03_docking/receptors_pdbqt/PqsR_4JVI.pdbqt",
        "04_scripts/python/step7_test_mmcif_to_pdbqt.py",
        f"99_logs_versions_checksums/run_logs/step7_LasR_2UV0_meeko_{DATE}.log",
        f"99_logs_versions_checksums/run_logs/step7_PqsR_4JVI_meeko_{DATE}.log",
        f"99_logs_versions_checksums/run_logs/step7_receptor_preparation_{DATE}.log",
    ]

    missing = [rel for rel in relpaths if not (root / rel).is_file()]
    if missing:
        fail("missing required files:\n" + "\n".join(missing))

    lasr_pdb = root / "00_validated_inputs/receptors/LasR_2UV0/LasR_2UV0_prepared_chainF.pdb"
    pqsr_pdb = root / "00_validated_inputs/receptors/PqsR_4JVI/PqsR_4JVI_prepared_dimer_AB.pdb"
    lasr_pdbqt = root / "03_docking/receptors_pdbqt/LasR_2UV0.pdbqt"
    pqsr_pdbqt = root / "03_docking/receptors_pdbqt/PqsR_4JVI.pdbqt"

    qc = {
        "lasr_pdb_atoms": count_records(lasr_pdb, ("ATOM  ",)),
        "pqsr_pdb_atoms": count_records(pqsr_pdb, ("ATOM  ",)),
        "lasr_pdbqt_atoms": count_records(lasr_pdbqt, ("ATOM", "HETATM")),
        "pqsr_pdbqt_atoms": count_records(pqsr_pdbqt, ("ATOM", "HETATM")),
        "lasr_charge": pdbqt_charge(lasr_pdbqt),
        "pqsr_charge": pdbqt_charge(pqsr_pdbqt),
    }
    expected = {
        "lasr_pdb_atoms": 1295,
        "pqsr_pdb_atoms": 3196,
        "lasr_pdbqt_atoms": 1565,
        "pqsr_pdbqt_atoms": 3910,
        "lasr_charge": -8.0,
        "pqsr_charge": -12.0,
    }
    for key, value in expected.items():
        if abs(qc[key] - value) > 1e-6:
            fail(f"QC mismatch for {key}: observed={qc[key]} expected={value}")

    versions_text = f"""Step 7 software versions
Date: {DATE}

Purpose:
Receptor preparation for LasR_2UV0 and PqsR_4JVI.
No docking was performed in Step 7.

Operating system:
{platform.platform()}

Python:
{sys.version.replace(chr(10), " ")}

Gemmi:
{pkg_version("gemmi")}

Meeko:
{pkg_version("meeko")}
Executable: /home/tomasz/miniforge3/envs/docking_env/bin/mk_prepare_receptor.py

RDKit:
{pkg_version("rdkit")}

Canonical preparation script:
04_scripts/python/step7_test_mmcif_to_pdbqt.py
SHA-256: {sha256(root / "04_scripts/python/step7_test_mmcif_to_pdbqt.py")}

Charge model:
Gasteiger charges generated by Meeko.

Receptor representation:
Rigid receptor PDBQT.
Polar hydrogens added by Meeko.
Unresolved polymer boundaries represented with Meeko --blunt_ends.
Waters excluded from the main receptor PDBQT files.
"""
    write_new(versions_path, versions_text)

    report_text = f"""# STEP 7 — Receptor preparation

**Date:** {DATE}  
**Status:** completed and QC-passed  
**Scope:** preparation of `LasR_2UV0` and `PqsR_4JVI` rigid receptor PDBQT files.  
**Docking:** not performed in this step.

## Canonical structural sources

The canonical archival inputs are the RCSB mmCIF files:

- `2UV0_rcsb_original.cif` — LasR ligand-binding domain, official RCSB resolution 1.8 Å.
- `4JVI_rcsb_original.cif` — PqsR/MvfR ligand-binding domain asymmetric unit, official RCSB resolution 2.9 Å.
- `4JVI_rcsb_biological_assembly1.cif` — biological dimer used for receptor preparation.

Official RCSB PDB exports were retained only as auxiliary compatibility records. The deterministic preparation script reads the canonical mmCIF files.

## LasR_2UV0 preparation

- Selected protomer: chain F, residues 5–168.
- Chain F was selected after structural QC because it had the most favorable pocket/protein crystallographic quality among the four protomers.
- Alternative conformations: altloc A retained for ARG137 and SER161; both are outside the orthosteric pocket.
- Four selenomethionines were converted by residue/atom identity only: MSE→MET and SE→SD at residues 26, 116, 144, and 153. Coordinates were not moved.
- OHN was removed from the docking receptor and preserved separately as `LasR_2UV0_OHN_chainF_reference.pdb`.
- OHN heavy-atom geometric centroid used as the prospective docking-box center: **(55.258, 27.138, 29.990) Å**.
- All waters were removed from the main receptor. The conserved bridging water HOH2134 was preserved separately for a possible sensitivity analysis; it is not part of the principal receptor ranking.
- LasR was prepared as a monomer because the OHN pocket is contained within one protomer; no partner-protomer residue lies within 6 Å of OHN.
- Meeko blunt ends: `F:5=0,F:168=2`.
- Prepared heavy-atom PDB: 1295 atoms.
- Final PDBQT: 1565 atoms, chain F, residues 5–168, net Gasteiger charge −8.0000.
- No OXT or artificial N-terminal H1/H2/H3/HT1/HT2/HT3 atoms were introduced.

## PqsR_4JVI preparation

- Biological assembly 1 was used because the QZN pocket includes the partner protomer; partner Leu153 approaches QZN to approximately 3.95 Å.
- mmCIF chains `A` and `A-2` were deterministically mapped to PDB chains A and B.
- Each protomer contains residues 94–296.
- QZN was removed from the docking receptor and preserved separately for both symmetry-related sites.
- QZN heavy-atom geometric centroids:
  - site A: **(−33.489, 56.865, 9.100) Å**;
  - site B: **(−32.502, 57.435, 29.370) Å**.
- Site A is designated as the primary reference site for later docking validation.
- All waters were removed from the main receptor. HOH501 was preserved separately for both symmetry-related sites but is excluded from the principal receptor.
- Four crystallographically unresolved backbone carbonyl oxygen atoms were reconstructed deterministically using the local CA(i)–C(i)–N(i+1) geometry and a C–O distance of 1.230 Å:
  - A:230 O = (−35.653, 57.066, −12.876) Å;
  - A:270 O = (−17.247, 65.535, 19.763) Å;
  - B:230 O = (−31.593, 59.410, 51.346) Å;
  - B:270 O = (−48.132, 47.705, 18.707) Å.
- These four oxygen atoms are modeled coordinates, not experimentally resolved atoms. No other protein heavy atom was moved.
- Meeko blunt ends: `A:94=0,A:296=2,B:94=0,B:296=2`.
- Prepared heavy-atom PDB: 3196 atoms.
- Final PDBQT: 3910 atoms, chains A/B, 203 residues per protomer, net Gasteiger charge −12.0000.
- No OXT or artificial N-terminal H1/H2/H3/HT1/HT2/HT3 atoms were introduced.

## Protonation and receptor representation

Meeko selected neutral histidine tautomers with protonation on NE2 for the receptor structures. No manual histidine-state override was applied because the audited histidine side chains are outside direct reference-ligand contact geometry. Receptors are rigid PDBQT structures with Gasteiger charges and Meeko-added polar hydrogens.

## Water policy

The principal receptor set is water-free for consistent Vina ranking across compounds. The highly conserved LasR bridging water was archived separately and may be evaluated only as a predefined sensitivity analysis. This water is not included in the main comparison. The weakly defined PqsR HOH501 was also excluded.

## QC summary

| Receptor | Prepared PDB atoms | PDBQT atoms | Chains | Residue range | Net charge |
|---|---:|---:|---|---|---:|
| LasR_2UV0 | 1295 | 1565 | F | 5–168 | −8.0000 |
| PqsR_4JVI | 3196 | 3910 | A, B | 94–296 each | −12.0000 |

Reference entities retained:

- OHN chain F: 21 atoms.
- LasR HOH2134: 1 atom.
- QZN site A: 22 atoms.
- QZN site B: 22 atoms.
- PqsR HOH501 site A: 1 atom.
- PqsR HOH501 site B: 1 atom.

## Reproducibility

Preparation script:

`04_scripts/python/step7_test_mmcif_to_pdbqt.py`

The script refuses to overwrite non-empty output directories, validates expected atom counts, chain identities, residue ranges, charges, reference-ligand centroids, MSE→MET conversion, modeled oxygen coordinates, and absence of terminal artifacts before reporting `STATUS: TEST PASSED`.

The SHA-256 manifest for all Step 7 inputs, outputs, scripts, logs, and documentation is:

`99_logs_versions_checksums/checksums/sha256_step7_receptors_{DATE}.txt`
"""
    write_new(report_path, report_text)

    manifest_relpaths = relpaths + [
        str(versions_path.relative_to(root)),
        str(report_path.relative_to(root)),
    ]
    manifest_lines = [
        f"{sha256(root / rel)}  {rel}"
        for rel in sorted(manifest_relpaths)
    ]
    write_new(manifest_path, "\n".join(manifest_lines) + "\n")

    print("=== STEP 7 DOCUMENTATION CREATED ===")
    for path in (versions_path, report_path, manifest_path):
        print(f"{path}\t{path.stat().st_size} bytes")
        print(f"SHA-256: {sha256(path)}")

    print("\n=== QC RECORDED ===")
    for key, value in qc.items():
        if isinstance(value, float):
            print(f"{key}={value:.4f}")
        else:
            print(f"{key}={value}")

    print("\n=== MANIFEST ENTRIES ===")
    print(len(manifest_lines))
    print("\nSTATUS: DOCUMENTATION PASSED")


if __name__ == "__main__":
    main()
