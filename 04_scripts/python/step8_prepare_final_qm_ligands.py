#!/usr/bin/env python3
"""
STEP 8C (ligand preparation) — dry-run/apply workflow.

Creates topology-preserving SDF files carrying the final QM coordinates and
Meeko/Gasteiger ligand PDBQT files for compound_08, compound_09, and compound_11.

Default mode is read-only dry-run. Use --apply to write repository outputs.

Core rules:
- final QM XYZ coordinates are transferred atom-for-atom onto audited MOL topology;
- no bond perception from Cartesian distances;
- no geometry optimization;
- neutral validated protomer/tautomer retained;
- explicit input hydrogens retained until Meeko's standard nonpolar-H merging;
- Gasteiger charges;
- conjugated bonds misclassified by default Meeko are rigidified by SMARTS:
  * compounds 08/09: hydrazone N–N in C=N–NH–C(=S)–NH2
  * compound 11: exocyclic NH2–thiadiazole C–N
- ring-N–scaffold sigma bond in compound_11 remains rotatable.

The script refuses to overwrite existing outputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import shutil
import subprocess
import tempfile
from collections import Counter
from datetime import date
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

ROOT = Path("/data/BIOFILM_INSILICO_SUPPORT")
RUN_DATE = "2026-07-16"

OUT_SDF_DIR = ROOT / "03_docking/ligands_sdf"
OUT_PDBQT_DIR = ROOT / "03_docking/ligands_pdbqt"
PREP_DIR = ROOT / "03_docking/ligand_preparation"
RUN_LOG = ROOT / f"99_logs_versions_checksums/run_logs/step8_ligand_pdbqt_preparation_{RUN_DATE}.md"
VERSION_FILE = ROOT / f"99_logs_versions_checksums/software_versions/step8_software_versions_{RUN_DATE}.txt"
CHECKSUM_FILE = ROOT / f"99_logs_versions_checksums/checksums/sha256_step8_ligands_{RUN_DATE}.txt"
PROTOCOL_FILE = PREP_DIR / "step8_ligand_preparation_protocol.json"

CASES = {
    "compound_08": {
        "formula": "C33H45N3O19S2",
        "atoms": 102,
        "formal_charge": 0,
        "xyz": ROOT / "01_qm_orca_xtb/compound_08/geometries/compound_08_wb97xd3_def2svp_cpcm_water_opt_min.xyz",
        "xyz_sha256": "ff0d4c0a6db6c794f345778133d1373355188c0614d9cbb808e77c7ed6e314fa",
        "topology": ROOT / "01_qm_orca_xtb/compound_08/geometries/compound_08_selected_for_xtb.mol",
        "smiles": ROOT / "00_validated_inputs/ligands/compound_08/compound_08.smiles",
        "rigid_smarts": "[C;X3](=[N;X2]-[N;X3]-[C;X3](=[S;X1])-[N;X3])",
        "rigid_bond_indices": ["2", "3"],
        "expected_torsdof": 20,
        "expected_pdbqt_atoms": 60,
        "rule": "Rigidify conjugated hydrazone N-N bond in C=N-NH-C(=S)-NH2.",
    },
    "compound_09": {
        "formula": "C18H25N3O9S2",
        "atoms": 57,
        "formal_charge": 0,
        "xyz": ROOT / "01_qm_orca_xtb/compound_09/geometries/compound_09_wb97xd3_def2svp_cpcm_water_minimum.xyz",
        "xyz_sha256": "bca8eb4e070e74bf0b81a48ff3bfe55406f289d9bcdeb2e4706dc5b61a0ee7aa",
        "topology": ROOT / "01_qm_orca_xtb/compound_09/geometries/compound_09_rdkit_ETKDGv3_best.mol",
        "smiles": ROOT / "00_validated_inputs/ligands/compound_09/compound_09.smiles",
        "rigid_smarts": "[C;X3](=[N;X2]-[N;X3]-[C;X3](=[S;X1])-[N;X3])",
        "rigid_bond_indices": ["2", "3"],
        "expected_torsdof": 8,
        "expected_pdbqt_atoms": 35,
        "rule": "Rigidify conjugated hydrazone N-N bond in C=N-NH-C(=S)-NH2.",
    },
    "compound_11": {
        "formula": "C8H9N3O3S2",
        "atoms": 25,
        "formal_charge": 0,
        "xyz": ROOT / "01_qm_orca_xtb/compound_11/geometries/compound_11_wb97xd3_def2svp_cpcm_water_minimum.xyz",
        "xyz_sha256": "8c5fa226efa04a7cbf24610c32e85ab9c0a549a68b1ec9bdcf5b27b94c4a273a",
        "topology": ROOT / "01_qm_orca_xtb/compound_11/geometries/compound_11_rdkit_ETKDGv3_best.mol",
        "smiles": ROOT / "00_validated_inputs/ligands/compound_11/compound_11.smiles",
        "rigid_smarts": "[N;X3;H2]-[c;R]",
        "rigid_bond_indices": ["1", "2"],
        "expected_torsdof": 1,
        "expected_pdbqt_atoms": 18,
        "rule": "Rigidify conjugated exocyclic NH2-thiadiazole C-N bond; retain ring-N-scaffold sigma rotor.",
    },
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_no_h_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)), canonical=True, isomericSmiles=True)


def read_xyz(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    declared = int(lines[0].strip())
    atoms = []
    for line in lines[2:]:
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 4:
            raise RuntimeError(f"Malformed XYZ line in {path}: {line}")
        atoms.append((parts[0], tuple(float(v) for v in parts[1:4])))
    if len(atoms) != declared:
        raise RuntimeError(f"{path}: declared {declared}, parsed {len(atoms)}")
    return atoms


def replace_coordinates(topology_mol: Chem.Mol, xyz_atoms):
    mol = Chem.Mol(topology_mol)
    if mol.GetNumAtoms() != len(xyz_atoms):
        raise RuntimeError("Topology/XYZ atom count mismatch")
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    for idx, (symbol, coords) in enumerate(xyz_atoms):
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() != symbol:
            raise RuntimeError(
                f"Element-order mismatch at atom {idx+1}: topology={atom.GetSymbol()} XYZ={symbol}"
            )
        conf.SetAtomPosition(idx, coords)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol


def max_coordinate_difference(mol: Chem.Mol, xyz_atoms) -> float:
    conf = mol.GetConformer()
    maximum = 0.0
    for idx, (_, coords) in enumerate(xyz_atoms):
        p = conf.GetAtomPosition(idx)
        diff = max(abs(p.x - coords[0]), abs(p.y - coords[1]), abs(p.z - coords[2]))
        maximum = max(maximum, diff)
    return maximum


def write_sdf(path: Path, mol: Chem.Mol, properties: dict[str, str]) -> None:
    out = Chem.Mol(mol)
    for key, value in properties.items():
        out.SetProp(key, str(value))
    writer = Chem.SDWriter(str(path))
    writer.SetForceV3000(True)
    writer.write(out)
    writer.close()


def read_single_sdf(path: Path) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
    mols = [m for m in supplier if m is not None]
    if len(mols) != 1:
        raise RuntimeError(f"{path}: expected one molecule, found {len(mols)}")
    return mols[0]


def run_command(cmd: list[str]) -> subprocess.CompletedProcess:
    proc = subprocess.run(cmd, text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            + " ".join(cmd)
            + f"\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return proc


def parse_torsdof(text: str) -> int:
    for line in text.splitlines():
        fields = line.split()
        if fields and fields[0] == "TORSDOF":
            return int(fields[1])
    raise RuntimeError("TORSDOF not found")


def parse_index_map(text: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    for line in text.splitlines():
        if line.startswith("REMARK INDEX MAP"):
            nums = [int(x) for x in re.findall(r"\d+", line[len("REMARK INDEX MAP"):])]
            if len(nums) % 2:
                raise RuntimeError(f"Odd INDEX MAP field count: {line}")
            for original_idx, pdbqt_serial in zip(nums[0::2], nums[1::2]):
                mapping[pdbqt_serial] = original_idx
    return mapping


def parse_pdbqt_atoms(text: str):
    atoms = {}
    charge_sum = 0.0
    bad_charge = False
    for line in text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            serial = int(line[6:11])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            fields = line.split()
            charge = float(fields[-2])
            if not (charge == charge and abs(charge) != float("inf")):
                bad_charge = True
            charge_sum += charge
            atoms[serial] = (x, y, z, charge, fields[-1])
    return atoms, charge_sum, bad_charge


def max_pdbqt_coordinate_difference(pdbqt_text: str, mol: Chem.Mol) -> float:
    mapping = parse_index_map(pdbqt_text)
    pdbqt_atoms, _, _ = parse_pdbqt_atoms(pdbqt_text)
    conf = mol.GetConformer()
    maximum = 0.0
    if not mapping:
        raise RuntimeError("PDBQT INDEX MAP missing")
    for serial, original_idx in mapping.items():
        if serial not in pdbqt_atoms:
            raise RuntimeError(f"PDBQT serial {serial} in INDEX MAP but no ATOM record")
        x, y, z, _, _ = pdbqt_atoms[serial]
        p = conf.GetAtomPosition(original_idx - 1)
        maximum = max(maximum, abs(x - p.x), abs(y - p.y), abs(z - p.z))
    return maximum


def package_version(name: str) -> str:
    try:
        from importlib.metadata import version
        return version(name)
    except Exception:
        return "not available"


def tool_output(args: list[str]) -> str:
    try:
        proc = subprocess.run(args, text=True, capture_output=True)
        return (proc.stdout + proc.stderr).strip()
    except Exception as exc:
        return f"ERROR: {exc}"


def output_paths():
    paths = []
    for name in CASES:
        paths += [
            OUT_SDF_DIR / f"{name}_final_qm_topology.sdf",
            OUT_PDBQT_DIR / f"{name}.pdbqt",
            PREP_DIR / f"{name}_meeko_prepare.log",
        ]
    paths += [PROTOCOL_FILE, RUN_LOG, VERSION_FILE, CHECKSUM_FILE]
    return paths


def validate_sources():
    errors = []
    summaries = {}
    for name, cfg in CASES.items():
        for key in ("xyz", "topology", "smiles"):
            if not cfg[key].is_file():
                errors.append(f"{name}: missing {key}: {cfg[key]}")
        if errors:
            continue

        actual_xyz_hash = sha256(cfg["xyz"])
        if actual_xyz_hash != cfg["xyz_sha256"]:
            errors.append(
                f"{name}: final XYZ hash mismatch: expected {cfg['xyz_sha256']} got {actual_xyz_hash}"
            )
            continue

        topology = Chem.MolFromMolFile(str(cfg["topology"]), removeHs=False, sanitize=True)
        validated = Chem.MolFromSmiles(cfg["smiles"].read_text(encoding="utf-8").strip())
        if topology is None or validated is None:
            errors.append(f"{name}: RDKit failed to read topology or validated SMILES")
            continue
        validated_h = Chem.AddHs(validated)
        xyz_atoms = read_xyz(cfg["xyz"])

        formula = rdMolDescriptors.CalcMolFormula(topology)
        charge = Chem.GetFormalCharge(topology)
        identity = canonical_no_h_smiles(topology) == canonical_no_h_smiles(validated_h)
        sequence = [a.GetSymbol() for a in topology.GetAtoms()] == [a[0] for a in xyz_atoms]

        if formula != cfg["formula"]:
            errors.append(f"{name}: formula mismatch {formula}")
        if topology.GetNumAtoms() != cfg["atoms"]:
            errors.append(f"{name}: atom-count mismatch {topology.GetNumAtoms()}")
        if charge != cfg["formal_charge"]:
            errors.append(f"{name}: formal-charge mismatch {charge}")
        if not identity:
            errors.append(f"{name}: stereochemical identity mismatch")
        if not sequence:
            errors.append(f"{name}: element-order mismatch")

        summaries[name] = {
            "formula": formula,
            "atoms": topology.GetNumAtoms(),
            "formal_charge": charge,
            "identity_match": identity,
            "element_sequence_match": sequence,
            "xyz_sha256": actual_xyz_hash,
            "topology_sha256": sha256(cfg["topology"]),
            "smiles_sha256": sha256(cfg["smiles"]),
        }

    if errors:
        raise RuntimeError("\n".join(errors))
    return summaries


def prepare_in_directory(workdir: Path, final_output: bool):
    results = {}
    for name, cfg in CASES.items():
        topology = Chem.MolFromMolFile(str(cfg["topology"]), removeHs=False, sanitize=True)
        validated = Chem.MolFromSmiles(cfg["smiles"].read_text(encoding="utf-8").strip())
        validated_h = Chem.AddHs(validated)
        xyz_atoms = read_xyz(cfg["xyz"])
        qm_mol = replace_coordinates(topology, xyz_atoms)

        sdf_path = workdir / f"{name}_final_qm_topology.sdf"
        pdbqt_path = workdir / f"{name}.pdbqt"
        log_path = workdir / f"{name}_meeko_prepare.log"

        props = {
            "compound_id": name,
            "source_qm_xyz": str(cfg["xyz"].relative_to(ROOT)),
            "source_qm_xyz_sha256": sha256(cfg["xyz"]),
            "source_topology_mol": str(cfg["topology"].relative_to(ROOT)),
            "source_topology_mol_sha256": sha256(cfg["topology"]),
            "validated_smiles_file": str(cfg["smiles"].relative_to(ROOT)),
            "validated_smiles_sha256": sha256(cfg["smiles"]),
            "geometry_level": "omegaB97X-D3/def2-SVP/CPCM(water)",
            "formal_charge": str(cfg["formal_charge"]),
            "protonation_tautomer_policy": "validated neutral input retained; no enumeration",
            "coordinate_transfer": "atom-for-atom from final QM XYZ; no distance-based bond perception",
            "meeko_charge_model": "gasteiger",
            "torsion_rule": cfg["rule"],
            "preparation_date": RUN_DATE,
        }
        write_sdf(sdf_path, qm_mol, props)
        reread = read_single_sdf(sdf_path)

        if rdMolDescriptors.CalcMolFormula(reread) != cfg["formula"]:
            raise RuntimeError(f"{name}: SDF formula failed after write/read")
        if Chem.GetFormalCharge(reread) != cfg["formal_charge"]:
            raise RuntimeError(f"{name}: SDF formal charge failed after write/read")
        if canonical_no_h_smiles(reread) != canonical_no_h_smiles(validated_h):
            raise RuntimeError(f"{name}: SDF stereochemical identity failed after write/read")
        sdf_coord_diff = max_coordinate_difference(reread, xyz_atoms)
        if sdf_coord_diff > 0.0005:
            raise RuntimeError(f"{name}: SDF coordinate precision loss too large: {sdf_coord_diff}")

        cmd = [
            "mk_prepare_ligand.py",
            "-i", str(sdf_path),
            "-o", str(pdbqt_path),
            "--charge_model", "gasteiger",
            "--add_index_map",
            "--rename_atoms",
            "-r", cfg["rigid_smarts"],
            "-b", *cfg["rigid_bond_indices"],
            "-v",
        ]
        proc = run_command(cmd)
        log_text = (
            "COMMAND:\n"
            + " ".join(cmd)
            + "\n\nSTDOUT:\n"
            + proc.stdout
            + "\nSTDERR:\n"
            + proc.stderr
        )
        log_path.write_text(log_text, encoding="utf-8")

        pdbqt_text = pdbqt_path.read_text(encoding="utf-8")
        torsdof = parse_torsdof(pdbqt_text)
        pdbqt_atoms, charge_sum, bad_charge = parse_pdbqt_atoms(pdbqt_text)
        pdbqt_coord_diff = max_pdbqt_coordinate_difference(pdbqt_text, reread)

        if torsdof != cfg["expected_torsdof"]:
            raise RuntimeError(f"{name}: TORSDOF expected {cfg['expected_torsdof']} got {torsdof}")
        if len(pdbqt_atoms) != cfg["expected_pdbqt_atoms"]:
            raise RuntimeError(
                f"{name}: PDBQT atoms expected {cfg['expected_pdbqt_atoms']} got {len(pdbqt_atoms)}"
            )
        if bad_charge:
            raise RuntimeError(f"{name}: NaN/Inf PDBQT charge")
        if abs(charge_sum) > 0.01:
            raise RuntimeError(f"{name}: PDBQT charge sum outside tolerance: {charge_sum}")
        if pdbqt_coord_diff > 0.0011:
            raise RuntimeError(
                f"{name}: PDBQT coordinate difference exceeds rounding tolerance: {pdbqt_coord_diff}"
            )

        results[name] = {
            "sdf": sdf_path,
            "pdbqt": pdbqt_path,
            "log": log_path,
            "formula": cfg["formula"],
            "input_atoms": cfg["atoms"],
            "pdbqt_atoms": len(pdbqt_atoms),
            "formal_charge": cfg["formal_charge"],
            "gasteiger_charge_sum": charge_sum,
            "torsdof": torsdof,
            "sdf_max_coordinate_difference_A": sdf_coord_diff,
            "pdbqt_max_coordinate_difference_A": pdbqt_coord_diff,
            "rigid_smarts": cfg["rigid_smarts"],
            "rigid_bond_indices": cfg["rigid_bond_indices"],
            "rule": cfg["rule"],
            "source_xyz": cfg["xyz"],
            "source_topology": cfg["topology"],
            "source_smiles": cfg["smiles"],
        }
    return results


def build_protocol(results):
    return {
        "step": 8,
        "date": RUN_DATE,
        "purpose": "Preparation of ligands for AutoDock Vina",
        "common_protocol": {
            "coordinate_source": "Final frozen STEP 3 omegaB97X-D3/def2-SVP/CPCM(water) XYZ geometry",
            "topology_source": "Audited RDKit MOL topology with identical atom order and validated stereochemical identity",
            "coordinate_transfer": "Atom-for-atom; no distance-based bond perception; no geometry optimization",
            "protonation_state": "Validated neutral state; total formal charge 0",
            "tautomer_handling": "No tautomer or protomer enumeration",
            "hydrogens": "Explicit QM hydrogens supplied; standard Meeko nonpolar-H merging used for PDBQT",
            "partial_charges": "Gasteiger generated by Meeko 0.7.1",
            "torsions": "Default Meeko torsions plus motif-based rigidification of conjugated exocyclic bonds misclassified as rotatable",
            "pdbqt_atom_mapping": "Meeko REMARK INDEX MAP retained",
            "pdbqt_atom_naming": "Meeko --rename_atoms",
        },
        "compound_rules": {
            name: {
                "source_qm_xyz": str(data["source_xyz"].relative_to(ROOT)),
                "source_topology": str(data["source_topology"].relative_to(ROOT)),
                "source_validated_smiles": str(data["source_smiles"].relative_to(ROOT)),
                "rigidify_smarts": data["rigid_smarts"],
                "rigidify_bond_indices_1_based_in_smarts": data["rigid_bond_indices"],
                "rule": data["rule"],
                "expected_torsdof": data["torsdof"],
            }
            for name, data in results.items()
        },
    }


def build_versions():
    lines = [
        "Step 8 software versions",
        f"Date: {RUN_DATE}",
        "",
        f"Operating system: {platform.platform()}",
        f"Python: {platform.python_version()}",
        f"Python executable: {shutil.which('python')}",
        f"RDKit: {Chem.rdBase.rdkitVersion}",
        f"Meeko: {package_version('meeko')}",
        f"Gemmi: {package_version('gemmi')}",
        f"mk_prepare_ligand.py: {shutil.which('mk_prepare_ligand.py')}",
        f"Open Babel: {tool_output(['obabel', '-V'])}",
        f"AutoDock Vina: {tool_output(['vina', '--version'])}",
        "",
        "Charge model: Gasteiger generated by Meeko.",
        "No docking was performed during ligand preparation.",
    ]
    return "\n".join(lines) + "\n"


def build_run_log(results):
    lines = [
        "# STEP 8 ligand PDBQT preparation",
        "",
        f"Date: {RUN_DATE}",
        "",
        "## Method",
        "",
        "- Final frozen STEP 3 QM coordinates were transferred atom-for-atom onto audited MOL topology.",
        "- No bond perception from distances and no geometry optimization were performed.",
        "- The validated neutral protomer/tautomer was retained for each compound.",
        "- Meeko assigned Gasteiger charges and merged standard nonpolar hydrogens.",
        "- REMARK INDEX MAP and deterministic atom renaming were retained.",
        "- Chemistry-aware SMARTS rules rigidified conjugated exocyclic bonds misclassified by default Meeko.",
        "",
        "## Results",
        "",
        "| Compound | Formula | Input atoms | PDBQT atoms | TORSDOF | Gasteiger charge sum | SDF max coordinate difference (Å) | PDBQT max coordinate difference (Å) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name, data in results.items():
        lines.append(
            f"| {name} | {data['formula']} | {data['input_atoms']} | {data['pdbqt_atoms']} | "
            f"{data['torsdof']} | {data['gasteiger_charge_sum']:.6f} | "
            f"{data['sdf_max_coordinate_difference_A']:.6f} | "
            f"{data['pdbqt_max_coordinate_difference_A']:.6f} |"
        )
    lines += [
        "",
        "## Torsion rules",
        "",
    ]
    for name, data in results.items():
        lines += [
            f"### {name}",
            "",
            f"- SMARTS: `{data['rigid_smarts']}`",
            f"- Bond indices: `{' '.join(data['rigid_bond_indices'])}`",
            f"- Rule: {data['rule']}",
            f"- Final TORSDOF: {data['torsdof']}",
            "",
        ]
    lines += [
        "## QC status",
        "",
        "All compounds passed formula, atom-order, stereochemical-identity, charge, TORSDOF, atom-count, and coordinate-preservation checks.",
        "",
        "Production docking was not performed in this step.",
        "",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true", help="write final repository outputs")
    args = parser.parse_args()

    source_summary = validate_sources()

    existing = [p for p in output_paths() if p.exists()]
    if existing:
        raise RuntimeError(
            "Refusing to overwrite existing outputs:\n" + "\n".join(str(p) for p in existing)
        )

    with tempfile.TemporaryDirectory(prefix="step8_ligprep_") as tmp:
        tmpdir = Path(tmp)
        dry_results = prepare_in_directory(tmpdir, final_output=False)

        print("SOURCE AUDIT PASSED")
        for name, summary in source_summary.items():
            print(
                f"{name}: formula={summary['formula']} atoms={summary['atoms']} "
                f"charge={summary['formal_charge']} identity={summary['identity_match']} "
                f"element_order={summary['element_sequence_match']}"
            )

        print("\nPREPARATION QC PASSED")
        for name, data in dry_results.items():
            print(
                f"{name}: PDBQT_ATOMS={data['pdbqt_atoms']} TORSDOF={data['torsdof']} "
                f"CHARGE_SUM={data['gasteiger_charge_sum']:.6f} "
                f"SDF_MAX_DXYZ={data['sdf_max_coordinate_difference_A']:.6f} A "
                f"PDBQT_MAX_DXYZ={data['pdbqt_max_coordinate_difference_A']:.6f} A"
            )

        if not args.apply:
            print("\nDRY-RUN PASSED")
            print("Repository was not modified.")
            print("Planned outputs:")
            for path in output_paths():
                print(f"  - {path}")
            return 0

    OUT_SDF_DIR.mkdir(parents=True, exist_ok=False)
    OUT_PDBQT_DIR.mkdir(parents=True, exist_ok=True)
    PREP_DIR.mkdir(parents=True, exist_ok=False)
    RUN_LOG.parent.mkdir(parents=True, exist_ok=True)
    VERSION_FILE.parent.mkdir(parents=True, exist_ok=True)
    CHECKSUM_FILE.parent.mkdir(parents=True, exist_ok=True)

    final_results = prepare_in_directory(PREP_DIR, final_output=True)

    # Move final SDF/PDBQT from staging preparation directory to their canonical destinations.
    for name, data in final_results.items():
        final_sdf = OUT_SDF_DIR / f"{name}_final_qm_topology.sdf"
        final_pdbqt = OUT_PDBQT_DIR / f"{name}.pdbqt"
        data["sdf"].replace(final_sdf)
        data["pdbqt"].replace(final_pdbqt)
        data["sdf"] = final_sdf
        data["pdbqt"] = final_pdbqt

    PROTOCOL_FILE.write_text(
        json.dumps(build_protocol(final_results), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    VERSION_FILE.write_text(build_versions(), encoding="utf-8")
    RUN_LOG.write_text(build_run_log(final_results), encoding="utf-8")

    checksum_targets = []
    for name, data in final_results.items():
        checksum_targets += [data["sdf"], data["pdbqt"], data["log"]]
    checksum_targets += [PROTOCOL_FILE, VERSION_FILE, RUN_LOG]

    lines = [f"{sha256(path)}  {path}" for path in checksum_targets]
    CHECKSUM_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\nAPPLY PASSED")
    for path in checksum_targets + [CHECKSUM_FILE]:
        print(f"{sha256(path)}  {path}")
    print("Backups have not yet been updated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
