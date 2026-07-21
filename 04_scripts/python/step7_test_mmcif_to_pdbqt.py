#!/usr/bin/env python3
"""
STEP 7 test utility: prepare LasR 2UV0 and PqsR 4JVI receptors from
canonical RCSB mmCIF files, convert to working PDB, run Meeko, and report QC.

Default behavior is non-destructive: all generated files are written to a new
temporary directory outside the project repository. Use --out-dir to preserve
outputs at a chosen location. This script does not perform docking.
"""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable

import gemmi


EXPECTED = {
    "lasr_pdb_atoms": 1295,
    "lasr_reference_atoms": 21,
    "lasr_water_atoms": 1,
    "lasr_pdbqt_atoms": 1565,
    "lasr_charge": -8.0,
    "pqsr_pdb_atoms": 3196,
    "pqsr_reference_atoms_per_chain": 22,
    "pqsr_water_atoms_per_chain": 1,
    "pqsr_pdbqt_atoms": 3910,
    "pqsr_charge": -12.0,
}


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def altloc_ok(atom: gemmi.Atom) -> bool:
    return str(atom.altloc) in ("", " ", "\x00", "A")


def atom_field(name: str, element: str) -> str:
    name = name.strip()
    element = element.strip().upper()
    if len(name) < 4 and len(element) == 1 and name[0].isalpha():
        return f" {name:<3s}"
    return f"{name:>4s}"


def pdb_atom_line(
    serial: int,
    record: str,
    atom_name: str,
    resname: str,
    chain: str,
    resnum: int,
    icode: str,
    xyz: tuple[float, float, float],
    occupancy: float,
    bfactor: float,
    element: str,
) -> str:
    ins = str(icode).strip()
    ins = ins[0] if ins else " "
    x, y, z = xyz
    return (
        f"{record:<6s}{serial:5d} {atom_field(atom_name, element)} "
        f"{resname:>3s} {chain:1s}{resnum:4d}{ins:1s}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}"
        f"          {element:>2s}  \n"
    )


def atom_data(atom: gemmi.Atom) -> dict:
    return {
        "name": atom.name.strip(),
        "xyz": (atom.pos.x, atom.pos.y, atom.pos.z),
        "occupancy": float(atom.occ),
        "bfactor": float(atom.b_iso),
        "element": atom.element.name,
    }


def unit(vector: tuple[float, float, float]) -> tuple[float, float, float]:
    norm = math.sqrt(sum(value * value for value in vector))
    if norm == 0:
        fail("zero vector encountered during carbonyl oxygen reconstruction")
    return tuple(value / norm for value in vector)


def centroid(coords: Iterable[tuple[float, float, float]]) -> tuple[float, float, float]:
    values = list(coords)
    if not values:
        fail("cannot calculate centroid of an empty coordinate set")
    return tuple(sum(point[i] for point in values) / len(values) for i in range(3))


def write_text_new(path: Path, content: str) -> None:
    if path.exists():
        fail(f"refusing to overwrite existing file: {path}")
    path.write_text(content, encoding="ascii")


def prepare_lasr(source_cif: Path, out_dir: Path) -> dict:
    structure = gemmi.read_structure(str(source_cif))
    if len(structure) != 1:
        fail(f"LasR: expected one model, found {len(structure)}")

    try:
        chain = structure[0]["F"]
    except KeyError:
        fail("LasR: chain F not found")

    protein_lines: list[str] = []
    ligand_lines: list[str] = []
    water_lines: list[str] = []
    ligand_xyz: list[tuple[float, float, float]] = []
    converted_mse: list[int] = []
    selected_altlocs: list[str] = []

    protein_serial = ligand_serial = water_serial = 1

    for residue in chain:
        resnum = residue.seqid.num
        resname = residue.name.strip()

        if residue.entity_type == gemmi.EntityType.Polymer and 5 <= resnum <= 168:
            output_resname = "MET" if resname == "MSE" else resname
            if resname == "MSE":
                converted_mse.append(resnum)

            for atom in residue:
                if not altloc_ok(atom):
                    continue
                data = atom_data(atom)
                if str(atom.altloc).strip():
                    selected_altlocs.append(f"F:{resnum}:{data['name']}:A")

                output_name = "SD" if resname == "MSE" and data["name"] == "SE" else data["name"]
                output_element = "S" if resname == "MSE" and data["name"] == "SE" else data["element"]

                protein_lines.append(
                    pdb_atom_line(
                        protein_serial,
                        "ATOM",
                        output_name,
                        output_resname,
                        "F",
                        resnum,
                        residue.seqid.icode,
                        data["xyz"],
                        data["occupancy"],
                        data["bfactor"],
                        output_element,
                    )
                )
                protein_serial += 1

        elif resname == "OHN":
            for atom in residue:
                if not altloc_ok(atom):
                    continue
                data = atom_data(atom)
                ligand_lines.append(
                    pdb_atom_line(
                        ligand_serial,
                        "HETATM",
                        data["name"],
                        "OHN",
                        "F",
                        resnum,
                        residue.seqid.icode,
                        data["xyz"],
                        data["occupancy"],
                        data["bfactor"],
                        data["element"],
                    )
                )
                ligand_xyz.append(data["xyz"])
                ligand_serial += 1

        elif resname == "HOH" and resnum == 2134:
            for atom in residue:
                data = atom_data(atom)
                water_lines.append(
                    pdb_atom_line(
                        water_serial,
                        "HETATM",
                        data["name"],
                        "HOH",
                        "F",
                        resnum,
                        residue.seqid.icode,
                        data["xyz"],
                        data["occupancy"],
                        data["bfactor"],
                        data["element"],
                    )
                )
                water_serial += 1

    write_text_new(out_dir / "LasR_2UV0_working.pdb", "".join(protein_lines) + "TER\nEND\n")
    write_text_new(
        out_dir / "LasR_2UV0_OHN_chainF_reference.pdb",
        "".join(ligand_lines) + "TER\nEND\n",
    )
    write_text_new(
        out_dir / "LasR_2UV0_HOH2134_chainF_reference.pdb",
        "".join(water_lines) + "TER\nEND\n",
    )

    unique_mse = sorted(set(converted_mse))
    if unique_mse != [26, 116, 144, 153]:
        fail(f"LasR: unexpected MSE residues converted: {unique_mse}")

    return {
        "protein_atoms": len(protein_lines),
        "ligand_atoms": len(ligand_lines),
        "water_atoms": len(water_lines),
        "ligand_centroid": centroid(ligand_xyz),
        "converted_mse": unique_mse,
        "selected_altlocs": selected_altlocs,
    }


def reconstruct_carbonyl_oxygen(
    atoms: list[dict],
    next_atoms: list[dict],
) -> tuple[float, float, float]:
    current = {atom["name"]: atom for atom in atoms}
    following = {atom["name"]: atom for atom in next_atoms}

    for name in ("CA", "C"):
        if name not in current:
            fail(f"missing required atom {name} for carbonyl reconstruction")
    if "N" not in following:
        fail("missing next-residue N for carbonyl reconstruction")

    ca = current["CA"]["xyz"]
    carbon = current["C"]["xyz"]
    next_n = following["N"]["xyz"]

    u_ca = unit(tuple(ca[i] - carbon[i] for i in range(3)))
    u_n = unit(tuple(next_n[i] - carbon[i] for i in range(3)))
    direction = unit(tuple(-(u_ca[i] + u_n[i]) for i in range(3)))
    return tuple(carbon[i] + 1.230 * direction[i] for i in range(3))


def prepare_pqsr(source_cif: Path, out_dir: Path) -> dict:
    structure = gemmi.read_structure(str(source_cif))
    if len(structure) != 1:
        fail(f"PqsR: expected one model, found {len(structure)}")

    model = structure[0]
    chain_map = {"A": "A", "A-2": "B"}
    protein_lines: list[str] = []
    ligand_lines: dict[str, list[str]] = {"A": [], "B": []}
    water_lines: dict[str, list[str]] = {"A": [], "B": []}
    ligand_xyz: dict[str, list[tuple[float, float, float]]] = {"A": [], "B": []}
    modeled_atoms: list[tuple[str, int, tuple[float, float, float]]] = []

    protein_serial = 1
    ligand_serial = {"A": 1, "B": 1}
    water_serial = {"A": 1, "B": 1}

    available_chains = {chain.name for chain in model}
    missing = set(chain_map) - available_chains
    if missing:
        fail(f"PqsR: missing assembly chains: {sorted(missing)}")

    for source_chain_name, output_chain in chain_map.items():
        chain = model[source_chain_name]
        residues = {
            residue.seqid.num: residue
            for residue in chain
            if residue.entity_type == gemmi.EntityType.Polymer
            and 94 <= residue.seqid.num <= 296
        }

        if sorted(residues) != list(range(94, 297)):
            fail(f"PqsR chain {output_chain}: residue range is incomplete")

        for resnum in range(94, 297):
            residue = residues[resnum]
            atoms = [atom_data(atom) for atom in residue if altloc_ok(atom)]
            atom_names = {atom["name"] for atom in atoms}

            if resnum in (230, 270):
                if "O" in atom_names:
                    fail(f"PqsR {output_chain}:{resnum} unexpectedly contains carbonyl O")
                next_atoms = [
                    atom_data(atom)
                    for atom in residues[resnum + 1]
                    if altloc_ok(atom)
                ]
                oxygen_xyz = reconstruct_carbonyl_oxygen(atoms, next_atoms)
                current = {atom["name"]: atom for atom in atoms}
                generated = {
                    "name": "O",
                    "xyz": oxygen_xyz,
                    "occupancy": current["C"]["occupancy"],
                    "bfactor": current["C"]["bfactor"],
                    "element": "O",
                }
                insert_index = next(i for i, atom in enumerate(atoms) if atom["name"] == "C") + 1
                atoms.insert(insert_index, generated)
                modeled_atoms.append((output_chain, resnum, oxygen_xyz))

            for atom in atoms:
                protein_lines.append(
                    pdb_atom_line(
                        protein_serial,
                        "ATOM",
                        atom["name"],
                        residue.name.strip(),
                        output_chain,
                        resnum,
                        residue.seqid.icode,
                        atom["xyz"],
                        atom["occupancy"],
                        atom["bfactor"],
                        atom["element"],
                    )
                )
                protein_serial += 1

        protein_lines.append("TER\n")

        for residue in chain:
            resname = residue.name.strip()
            resnum = residue.seqid.num

            if resname == "QZN":
                for atom in residue:
                    if not altloc_ok(atom):
                        continue
                    data = atom_data(atom)
                    ligand_lines[output_chain].append(
                        pdb_atom_line(
                            ligand_serial[output_chain],
                            "HETATM",
                            data["name"],
                            "QZN",
                            output_chain,
                            resnum,
                            residue.seqid.icode,
                            data["xyz"],
                            data["occupancy"],
                            data["bfactor"],
                            data["element"],
                        )
                    )
                    ligand_xyz[output_chain].append(data["xyz"])
                    ligand_serial[output_chain] += 1

            elif resname == "HOH" and resnum == 501:
                for atom in residue:
                    data = atom_data(atom)
                    water_lines[output_chain].append(
                        pdb_atom_line(
                            water_serial[output_chain],
                            "HETATM",
                            data["name"],
                            "HOH",
                            output_chain,
                            resnum,
                            residue.seqid.icode,
                            data["xyz"],
                            data["occupancy"],
                            data["bfactor"],
                            data["element"],
                        )
                    )
                    water_serial[output_chain] += 1

    write_text_new(out_dir / "PqsR_4JVI_working.pdb", "".join(protein_lines) + "END\n")

    for chain in ("A", "B"):
        write_text_new(
            out_dir / f"PqsR_4JVI_QZN_chain{chain}_reference.pdb",
            "".join(ligand_lines[chain]) + "TER\nEND\n",
        )
        write_text_new(
            out_dir / f"PqsR_4JVI_HOH501_chain{chain}_reference.pdb",
            "".join(water_lines[chain]) + "TER\nEND\n",
        )

    return {
        "protein_atoms": sum(line.startswith("ATOM") for line in protein_lines),
        "ligand_atoms": {chain: len(ligand_lines[chain]) for chain in ("A", "B")},
        "water_atoms": {chain: len(water_lines[chain]) for chain in ("A", "B")},
        "ligand_centroids": {
            chain: centroid(ligand_xyz[chain]) for chain in ("A", "B")
        },
        "modeled_atoms": modeled_atoms,
    }


def run_meeko(
    meeko_executable: Path,
    input_pdb: Path,
    output_pdbqt: Path,
    blunt_ends: str,
    log_path: Path,
) -> None:
    command = [
        str(meeko_executable),
        "--read_pdb",
        str(input_pdb),
        "-p",
        str(output_pdbqt),
        "--charge_model",
        "gasteiger",
        "--blunt_ends",
        blunt_ends,
    ]
    result = subprocess.run(command, text=True, capture_output=True)
    log_path.write_text(
        "COMMAND:\n"
        + " ".join(command)
        + "\n\nSTDOUT:\n"
        + result.stdout
        + "\nSTDERR:\n"
        + result.stderr,
        encoding="utf-8",
    )
    if result.returncode != 0:
        fail(f"Meeko failed for {input_pdb.name}; see {log_path}")


def read_pdbqt(path: Path) -> list[dict]:
    atoms: list[dict] = []
    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            fields = line.split()
            atoms.append(
                {
                    "name": line[12:16].strip(),
                    "resname": line[17:20].strip(),
                    "chain": line[21].strip(),
                    "resnum": int(line[22:26]),
                    "charge": float(fields[-2]),
                    "type": fields[-1],
                }
            )
    return atoms


def qc_pdbqt(
    path: Path,
    expected_atoms: int,
    expected_chains: set[str],
    expected_charge: float,
    expected_ranges: dict[str, tuple[int, int]],
) -> dict:
    atoms = read_pdbqt(path)
    if len(atoms) != expected_atoms:
        fail(f"{path.name}: atom count {len(atoms)} != expected {expected_atoms}")

    chains = {atom["chain"] for atom in atoms}
    if chains != expected_chains:
        fail(f"{path.name}: chains {sorted(chains)} != expected {sorted(expected_chains)}")

    charge = sum(atom["charge"] for atom in atoms)
    if abs(charge - expected_charge) > 0.01:
        fail(f"{path.name}: charge {charge:.4f} != expected {expected_charge:.4f}")

    forbidden = {"OXT", "H1", "H2", "H3", "HT1", "HT2", "HT3"}
    found_forbidden = sorted({atom["name"] for atom in atoms if atom["name"] in forbidden})
    if found_forbidden:
        fail(f"{path.name}: forbidden terminal atoms present: {found_forbidden}")

    ranges = {}
    for chain, (start, end) in expected_ranges.items():
        observed = sorted({atom["resnum"] for atom in atoms if atom["chain"] == chain})
        expected = list(range(start, end + 1))
        if observed != expected:
            fail(f"{path.name}: incomplete residue range for chain {chain}")
        ranges[chain] = (observed[0], observed[-1], len(observed))

    return {
        "atoms": len(atoms),
        "chains": sorted(chains),
        "charge": charge,
        "ranges": ranges,
        "forbidden_terminal_atoms": found_forbidden,
    }


def verify_pre_meeko(lasr: dict, pqsr: dict) -> None:
    checks = [
        (lasr["protein_atoms"], EXPECTED["lasr_pdb_atoms"], "LasR working PDB atoms"),
        (lasr["ligand_atoms"], EXPECTED["lasr_reference_atoms"], "LasR OHN atoms"),
        (lasr["water_atoms"], EXPECTED["lasr_water_atoms"], "LasR reference water atoms"),
        (pqsr["protein_atoms"], EXPECTED["pqsr_pdb_atoms"], "PqsR working PDB atoms"),
    ]
    for observed, expected, label in checks:
        if observed != expected:
            fail(f"{label}: {observed} != expected {expected}")

    for chain in ("A", "B"):
        if pqsr["ligand_atoms"][chain] != EXPECTED["pqsr_reference_atoms_per_chain"]:
            fail(f"PqsR QZN chain {chain}: unexpected atom count")
        if pqsr["water_atoms"][chain] != EXPECTED["pqsr_water_atoms_per_chain"]:
            fail(f"PqsR HOH501 chain {chain}: unexpected atom count")

    modeled_keys = {(chain, resnum) for chain, resnum, _ in pqsr["modeled_atoms"]}
    expected_keys = {("A", 230), ("A", 270), ("B", 230), ("B", 270)}
    if modeled_keys != expected_keys:
        fail(f"PqsR modeled atoms mismatch: {sorted(modeled_keys)}")


def print_summary(out_dir: Path, lasr: dict, pqsr: dict, lasr_qc: dict, pqsr_qc: dict) -> None:
    print("\n=== STEP 7 TEST: mmCIF -> working PDB -> PDBQT ===")
    print(f"Output directory: {out_dir}")

    print("\nLasR_2UV0")
    print(f"  working PDB atoms: {lasr['protein_atoms']}")
    print(f"  OHN reference atoms: {lasr['ligand_atoms']}")
    print(f"  conserved water atoms: {lasr['water_atoms']}")
    print("  OHN centroid: {:.3f} {:.3f} {:.3f}".format(*lasr["ligand_centroid"]))
    print(f"  MSE->MET residues: {lasr['converted_mse']}")
    print(f"  selected altloc A atoms: {len(lasr['selected_altlocs'])}")
    print(f"  PDBQT atoms: {lasr_qc['atoms']}")
    print(f"  PDBQT charge: {lasr_qc['charge']:.4f}")
    print(f"  PDBQT chains: {lasr_qc['chains']}")
    print(f"  PDBQT residue ranges: {lasr_qc['ranges']}")

    print("\nPqsR_4JVI")
    print(f"  working PDB atoms: {pqsr['protein_atoms']}")
    print(f"  QZN atoms per chain: {pqsr['ligand_atoms']}")
    print(f"  HOH501 atoms per chain: {pqsr['water_atoms']}")
    for chain in ("A", "B"):
        print(
            "  QZN chain {} centroid: {:.3f} {:.3f} {:.3f}".format(
                chain, *pqsr["ligand_centroids"][chain]
            )
        )
    print("  modeled carbonyl O atoms:")
    for chain, resnum, xyz in pqsr["modeled_atoms"]:
        print(f"    {chain}:{resnum} O = {xyz[0]:.3f} {xyz[1]:.3f} {xyz[2]:.3f}")
    print(f"  PDBQT atoms: {pqsr_qc['atoms']}")
    print(f"  PDBQT charge: {pqsr_qc['charge']:.4f}")
    print(f"  PDBQT chains: {pqsr_qc['chains']}")
    print(f"  PDBQT residue ranges: {pqsr_qc['ranges']}")

    print("\nSTATUS: TEST PASSED")
    print("No docking was performed.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        default="/data/BIOFILM_INSILICO_SUPPORT",
        help="BIOFILM_INSILICO_SUPPORT root directory",
    )
    parser.add_argument(
        "--meeko",
        default="/home/tomasz/miniforge3/envs/docking_env/bin/mk_prepare_receptor.py",
        help="Path to mk_prepare_receptor.py",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="New empty output directory. Default: create a new directory under /tmp.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = Path(args.root).resolve()
    meeko = Path(args.meeko).resolve()

    if not root.is_dir():
        fail(f"project root not found: {root}")
    if not meeko.is_file():
        fail(f"Meeko executable not found: {meeko}")

    lasr_cif = root / "00_validated_inputs/receptors/LasR_2UV0/2UV0_rcsb_original.cif"
    pqsr_cif = (
        root
        / "00_validated_inputs/receptors/PqsR_4JVI/4JVI_rcsb_biological_assembly1.cif"
    )
    for source in (lasr_cif, pqsr_cif):
        if not source.is_file():
            fail(f"source mmCIF not found: {source}")

    if args.out_dir:
        out_dir = Path(args.out_dir).resolve()
        if out_dir.exists():
            if any(out_dir.iterdir()):
                fail(f"output directory is not empty: {out_dir}")
        else:
            out_dir.mkdir(parents=True)
    else:
        out_dir = Path(tempfile.mkdtemp(prefix="step7_receptor_test_"))

    lasr = prepare_lasr(lasr_cif, out_dir)
    pqsr = prepare_pqsr(pqsr_cif, out_dir)
    verify_pre_meeko(lasr, pqsr)

    lasr_pdbqt = out_dir / "LasR_2UV0.pdbqt"
    pqsr_pdbqt = out_dir / "PqsR_4JVI.pdbqt"

    run_meeko(
        meeko,
        out_dir / "LasR_2UV0_working.pdb",
        lasr_pdbqt,
        "F:5=0,F:168=2",
        out_dir / "LasR_2UV0_meeko.log",
    )
    run_meeko(
        meeko,
        out_dir / "PqsR_4JVI_working.pdb",
        pqsr_pdbqt,
        "A:94=0,A:296=2,B:94=0,B:296=2",
        out_dir / "PqsR_4JVI_meeko.log",
    )

    lasr_qc = qc_pdbqt(
        lasr_pdbqt,
        EXPECTED["lasr_pdbqt_atoms"],
        {"F"},
        EXPECTED["lasr_charge"],
        {"F": (5, 168)},
    )
    pqsr_qc = qc_pdbqt(
        pqsr_pdbqt,
        EXPECTED["pqsr_pdbqt_atoms"],
        {"A", "B"},
        EXPECTED["pqsr_charge"],
        {"A": (94, 296), "B": (94, 296)},
    )

    print_summary(out_dir, lasr, pqsr, lasr_qc, pqsr_qc)


if __name__ == "__main__":
    main()
