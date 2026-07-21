#!/usr/bin/env python3
"""
STEP 8D — finalize validated 30 × 30 × 30 Å docking boxes.

Default mode is a read-only dry-run. Use --apply only after the dry-run passes.

The script:
- recomputes crystallographic reference-ligand centroids;
- confirms all reference heavy atoms are inside the box;
- confirms all receptor heavy atoms within 6 Å of the reference ligand are inside;
- confirms compounds 8, 9, 11 and the reference ligand fit at the box center
  under arbitrary rotation using their maximum centroid radius;
- writes only box-definition files and documentation;
- does not run Vina and does not create docking poses.

The same receptor-specific box is intended for reference redocking and later
compound docking.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import shutil
from pathlib import Path

import numpy as np
from rdkit import Chem

ROOT = Path("/data/BIOFILM_INSILICO_SUPPORT")
RUN_DATE = "2026-07-16"
BOX_SIZE = np.array([30.0, 30.0, 30.0], dtype=float)
HALF_BOX = BOX_SIZE / 2.0
POCKET_CUTOFF_A = 6.0
CENTER_TOLERANCE_A = 0.001

BOX_DIR = ROOT / "03_docking/configs/box_definitions"
DOC_DIR = ROOT / "03_docking/box_definition"
SCRIPT_DEST = ROOT / "04_scripts/python/step8d_define_docking_boxes.py"

PROTOCOL_FILE = DOC_DIR / "step8d_docking_box_protocol.json"
RUN_LOG = ROOT / f"99_logs_versions_checksums/run_logs/step8d_docking_box_definition_{RUN_DATE}.md"
VERSION_FILE = ROOT / f"99_logs_versions_checksums/software_versions/step8d_software_versions_{RUN_DATE}.txt"
CHECKSUM_FILE = ROOT / f"99_logs_versions_checksums/checksums/sha256_step8d_docking_boxes_{RUN_DATE}.txt"

LIGANDS = {
    "compound_08": ROOT / "03_docking/ligands_sdf/compound_08_final_qm_topology.sdf",
    "compound_09": ROOT / "03_docking/ligands_sdf/compound_09_final_qm_topology.sdf",
    "compound_11": ROOT / "03_docking/ligands_sdf/compound_11_final_qm_topology.sdf",
}

SYSTEMS = {
    "LasR_2UV0": {
        "receptor_pdb": ROOT / "00_validated_inputs/receptors/LasR_2UV0/LasR_2UV0_prepared_chainF.pdb",
        "receptor_pdbqt": ROOT / "03_docking/receptors_pdbqt/LasR_2UV0.pdbqt",
        "reference_pdb": ROOT / "00_validated_inputs/receptors/LasR_2UV0/LasR_2UV0_OHN_chainF_reference.pdb",
        "reference_sdf": ROOT / "03_docking/reference_ligands_sdf/LasR_2UV0_OHN_chainF_reference_prepared.sdf",
        "reference_pdbqt": ROOT / "03_docking/reference_ligands_pdbqt/LasR_2UV0_OHN_chainF_reference.pdbqt",
        "reference_name": "OHN",
        "expected_center": np.array([55.258, 27.138, 29.990], dtype=float),
        "box_file": BOX_DIR / "LasR_2UV0_box_30A.txt",
    },
    "PqsR_4JVI": {
        "receptor_pdb": ROOT / "00_validated_inputs/receptors/PqsR_4JVI/PqsR_4JVI_prepared_dimer_AB.pdb",
        "receptor_pdbqt": ROOT / "03_docking/receptors_pdbqt/PqsR_4JVI.pdbqt",
        "reference_pdb": ROOT / "00_validated_inputs/receptors/PqsR_4JVI/PqsR_4JVI_QZN_chainA_reference.pdb",
        "reference_sdf": ROOT / "03_docking/reference_ligands_sdf/PqsR_4JVI_QZN_siteA_reference_prepared.sdf",
        "reference_pdbqt": ROOT / "03_docking/reference_ligands_pdbqt/PqsR_4JVI_QZN_siteA_reference.pdbqt",
        "reference_name": "QZN_siteA",
        "expected_center": np.array([-33.489, 56.865, 9.100], dtype=float),
        "box_file": BOX_DIR / "PqsR_4JVI_box_30A.txt",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pdb_heavy_coordinates(path: Path, record_types: tuple[str, ...]) -> np.ndarray:
    coordinates = []

    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(record_types):
            continue

        element = line[76:78].strip().upper()
        if not element:
            atom_name = "".join(character for character in line[12:16] if character.isalpha())
            element = atom_name[:1].upper()

        if element == "H":
            continue

        coordinates.append(
            [
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ]
        )

    if not coordinates:
        raise RuntimeError(f"No heavy atoms found in {path}")

    return np.asarray(coordinates, dtype=float)


def sdf_heavy_coordinates(path: Path) -> np.ndarray:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
    molecules = [molecule for molecule in supplier if molecule is not None]

    if len(molecules) != 1:
        raise RuntimeError(f"{path}: expected one molecule, found {len(molecules)}")

    molecule = molecules[0]
    conformer = molecule.GetConformer()

    coordinates = [
        list(conformer.GetAtomPosition(atom.GetIdx()))
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() > 1
    ]

    if not coordinates:
        raise RuntimeError(f"No heavy atoms found in {path}")

    return np.asarray(coordinates, dtype=float)


def centroid(coordinates: np.ndarray) -> np.ndarray:
    return coordinates.mean(axis=0)


def maximum_centroid_radius(coordinates: np.ndarray) -> float:
    center = centroid(coordinates)
    return float(np.linalg.norm(coordinates - center, axis=1).max())


def inside_box(coordinates: np.ndarray, center: np.ndarray) -> np.ndarray:
    return np.all(np.abs(coordinates - center) <= HALF_BOX + 1e-9, axis=1)


def minimum_face_clearance(coordinates: np.ndarray, center: np.ndarray) -> float:
    return float(np.min(HALF_BOX - np.abs(coordinates - center)))


def format_vector(values: np.ndarray) -> str:
    return " ".join(f"{float(value):.3f}" for value in values)


def audit_system(system_name: str, configuration: dict, ligand_radii: dict[str, float]) -> dict:
    reference_coordinates = pdb_heavy_coordinates(
        configuration["reference_pdb"],
        ("ATOM  ", "HETATM"),
    )
    center = centroid(reference_coordinates)
    center_delta = np.abs(center - configuration["expected_center"])

    receptor_coordinates = pdb_heavy_coordinates(
        configuration["receptor_pdb"],
        ("ATOM  ",),
    )
    distances = np.linalg.norm(
        receptor_coordinates[:, None, :] - reference_coordinates[None, :, :],
        axis=2,
    )
    pocket_coordinates = receptor_coordinates[
        np.min(distances, axis=1) <= POCKET_CUTOFF_A
    ]

    reference_radius = maximum_centroid_radius(
        sdf_heavy_coordinates(configuration["reference_sdf"])
    )
    all_radii = {
        **ligand_radii,
        configuration["reference_name"]: reference_radius,
    }
    largest_ligand = max(all_radii, key=all_radii.get)
    largest_radius = all_radii[largest_ligand]

    checks = {
        "center_matches_audited_value": bool(np.max(center_delta) <= CENTER_TOLERANCE_A),
        "reference_inside_box": bool(np.all(inside_box(reference_coordinates, center))),
        "pocket_inside_box": bool(np.all(inside_box(pocket_coordinates, center))),
        "all_ligands_rotation_safe_at_center": bool(largest_radius <= float(np.min(HALF_BOX))),
        "box_volume_is_27000_A3": math.isclose(float(np.prod(BOX_SIZE)), 27000.0),
    }

    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise RuntimeError(f"{system_name}: failed checks: {failed}")

    return {
        "system": system_name,
        "center_A": [round(float(value), 3) for value in center],
        "size_A": [30.0, 30.0, 30.0],
        "volume_A3": 27000.0,
        "reference_name": configuration["reference_name"],
        "reference_heavy_atoms": int(len(reference_coordinates)),
        "pocket_cutoff_A": POCKET_CUTOFF_A,
        "pocket_heavy_atoms": int(len(pocket_coordinates)),
        "reference_minimum_face_clearance_A": minimum_face_clearance(
            reference_coordinates,
            center,
        ),
        "pocket_minimum_face_clearance_A": minimum_face_clearance(
            pocket_coordinates,
            center,
        ),
        "ligand_centroid_radii_A": all_radii,
        "largest_ligand_centroid_radius_A": largest_radius,
        "largest_ligand": largest_ligand,
        "rotation_clearance_at_center_A": float(np.min(HALF_BOX)) - largest_radius,
        "checks": checks,
    }


def box_text(result: dict) -> str:
    center_x, center_y, center_z = result["center_A"]
    size_x, size_y, size_z = result["size_A"]

    return "\n".join(
        [
            f"center_x = {center_x:.3f}",
            f"center_y = {center_y:.3f}",
            f"center_z = {center_z:.3f}",
            f"size_x = {size_x:.3f}",
            f"size_y = {size_y:.3f}",
            f"size_z = {size_z:.3f}",
            "",
        ]
    )


def protocol(results: dict[str, dict]) -> dict:
    systems = {}

    for system_name, result in results.items():
        configuration = SYSTEMS[system_name]
        systems[system_name] = {
            **result,
            "box_file": str(configuration["box_file"].relative_to(ROOT)),
            "receptor_pdb": str(configuration["receptor_pdb"].relative_to(ROOT)),
            "receptor_pdb_sha256": sha256(configuration["receptor_pdb"]),
            "receptor_pdbqt": str(configuration["receptor_pdbqt"].relative_to(ROOT)),
            "receptor_pdbqt_sha256": sha256(configuration["receptor_pdbqt"]),
            "reference_pdb": str(configuration["reference_pdb"].relative_to(ROOT)),
            "reference_pdb_sha256": sha256(configuration["reference_pdb"]),
            "reference_sdf": str(configuration["reference_sdf"].relative_to(ROOT)),
            "reference_sdf_sha256": sha256(configuration["reference_sdf"]),
            "reference_pdbqt": str(configuration["reference_pdbqt"].relative_to(ROOT)),
            "reference_pdbqt_sha256": sha256(configuration["reference_pdbqt"]),
        }

    return {
        "step": "8D",
        "date": RUN_DATE,
        "purpose": "Validated receptor-specific docking-box definitions for crystallographic redocking and later docking of compounds 8, 9, and 11",
        "box_policy": {
            "shape": "cube",
            "size_A": [30.0, 30.0, 30.0],
            "volume_A3": 27000.0,
            "center_definition": "Heavy-atom geometric centroid of the crystallographic reference ligand",
            "pocket_validation": "All receptor heavy atoms within 6.0 Å of any reference-ligand heavy atom must lie inside the box",
            "ligand_validation": "The maximum centroid radius among compounds 8, 9, 11 and the receptor-specific reference ligand must fit at the box center under arbitrary rotation",
            "common_use": "The same receptor-specific box must be used for reference redocking and subsequent compound docking",
            "vina_parameters": "Not defined in STEP 8D",
            "docking_status": "No redocking or production docking performed in STEP 8D",
        },
        "systems": systems,
    }


def run_log(results: dict[str, dict]) -> str:
    lines = [
        "# STEP 8D docking-box definition",
        "",
        f"Date: {RUN_DATE}",
        "",
        "## Final box policy",
        "",
        "- Box shape: cube.",
        "- Box dimensions: 30.000 × 30.000 × 30.000 Å.",
        "- Box volume: 27,000 Å³.",
        "- Center: crystallographic reference-ligand heavy-atom geometric centroid.",
        "- The same receptor-specific box will be used for reference redocking and later docking of compounds 8, 9, and 11.",
        "- All receptor heavy atoms within 6.0 Å of the crystallographic reference ligand were required to lie inside the box.",
        "- All ligands were required to fit at the box center under arbitrary rotation, using maximum centroid radius.",
        "- Vina search parameters were not defined and docking was not run.",
        "",
        "## Validated systems",
        "",
        "| System | Center x | Center y | Center z | Size x | Size y | Size z | Pocket atoms | Pocket clearance | Largest ligand | Radius | Rotation clearance |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|",
    ]

    for system_name, result in results.items():
        center = result["center_A"]
        size = result["size_A"]
        lines.append(
            f"| {system_name} | {center[0]:.3f} | {center[1]:.3f} | "
            f"{center[2]:.3f} | {size[0]:.3f} | {size[1]:.3f} | "
            f"{size[2]:.3f} | {result['pocket_heavy_atoms']} | "
            f"{result['pocket_minimum_face_clearance_A']:.3f} Å | "
            f"{result['largest_ligand']} | "
            f"{result['largest_ligand_centroid_radius_A']:.3f} Å | "
            f"{result['rotation_clearance_at_center_A']:.3f} Å |"
        )

    lines += [
        "",
        "## QC status",
        "",
        "Both receptor-specific boxes passed center, reference-ligand containment, pocket containment, arbitrary-rotation ligand-fit, and volume checks.",
        "",
        "No redocking or production docking was performed.",
        "",
    ]

    return "\n".join(lines)


def versions_text() -> str:
    return "\n".join(
        [
            "STEP 8D software versions",
            f"Date: {RUN_DATE}",
            "",
            f"Operating system: {platform.platform()}",
            f"Python: {platform.python_version()}",
            f"RDKit: {Chem.rdBase.rdkitVersion}",
            f"NumPy: {np.__version__}",
            "",
            "No Vina execution occurred in STEP 8D.",
            "",
        ]
    )


def planned_outputs() -> list[Path]:
    return [
        SYSTEMS["LasR_2UV0"]["box_file"],
        SYSTEMS["PqsR_4JVI"]["box_file"],
        PROTOCOL_FILE,
        RUN_LOG,
        VERSION_FILE,
        SCRIPT_DEST,
        CHECKSUM_FILE,
    ]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--apply",
        action="store_true",
        help="write validated STEP 8D box definitions and documentation",
    )
    arguments = parser.parse_args()

    source_script = Path(__file__).resolve()
    if not source_script.is_file():
        raise RuntimeError("Cannot resolve the running script path")

    required_paths = list(LIGANDS.values())
    for configuration in SYSTEMS.values():
        required_paths.extend(
            [
                configuration["receptor_pdb"],
                configuration["receptor_pdbqt"],
                configuration["reference_pdb"],
                configuration["reference_sdf"],
                configuration["reference_pdbqt"],
            ]
        )

    for path in required_paths:
        if not path.is_file():
            raise FileNotFoundError(path)

    existing = [path for path in planned_outputs() if path.exists()]
    if existing:
        raise RuntimeError(
            "Refusing to overwrite existing STEP 8D outputs:\n"
            + "\n".join(str(path) for path in existing)
        )

    ligand_radii = {
        name: maximum_centroid_radius(sdf_heavy_coordinates(path))
        for name, path in LIGANDS.items()
    }

    results = {
        system_name: audit_system(system_name, configuration, ligand_radii)
        for system_name, configuration in SYSTEMS.items()
    }

    print("STEP 8D BOX AUDIT PASSED")
    for system_name, result in results.items():
        print(
            f"{system_name}: "
            f"CENTER={format_vector(np.asarray(result['center_A']))} A "
            f"SIZE={format_vector(BOX_SIZE)} A "
            f"POCKET_ATOMS={result['pocket_heavy_atoms']} "
            f"POCKET_CLEARANCE={result['pocket_minimum_face_clearance_A']:.3f} A "
            f"LARGEST_RADIUS={result['largest_ligand_centroid_radius_A']:.3f} A "
            f"ROTATION_CLEARANCE={result['rotation_clearance_at_center_A']:.3f} A"
        )

    if not arguments.apply:
        print("\nDRY-RUN PASSED")
        print("Repository was not modified.")
        print("Planned outputs:")
        for path in planned_outputs():
            print(f"  - {path}")
        return 0

    BOX_DIR.mkdir(parents=True, exist_ok=True)
    DOC_DIR.mkdir(parents=True, exist_ok=True)
    RUN_LOG.parent.mkdir(parents=True, exist_ok=True)
    VERSION_FILE.parent.mkdir(parents=True, exist_ok=True)
    CHECKSUM_FILE.parent.mkdir(parents=True, exist_ok=True)
    SCRIPT_DEST.parent.mkdir(parents=True, exist_ok=True)

    for system_name, result in results.items():
        SYSTEMS[system_name]["box_file"].write_text(
            box_text(result),
            encoding="utf-8",
        )

    PROTOCOL_FILE.write_text(
        json.dumps(protocol(results), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    RUN_LOG.write_text(run_log(results), encoding="utf-8")
    VERSION_FILE.write_text(versions_text(), encoding="utf-8")
    shutil.copy2(source_script, SCRIPT_DEST)

    checksum_targets = [
        SYSTEMS["LasR_2UV0"]["box_file"],
        SYSTEMS["PqsR_4JVI"]["box_file"],
        PROTOCOL_FILE,
        RUN_LOG,
        VERSION_FILE,
        SCRIPT_DEST,
    ]

    CHECKSUM_FILE.write_text(
        "\n".join(
            f"{sha256(path)}  {path}"
            for path in checksum_targets
        )
        + "\n",
        encoding="utf-8",
    )

    print("\nAPPLY PASSED")
    for path in checksum_targets + [CHECKSUM_FILE]:
        print(f"{sha256(path)}  {path}")

    print("Backups have not yet been updated.")
    print("No redocking or production docking was performed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
