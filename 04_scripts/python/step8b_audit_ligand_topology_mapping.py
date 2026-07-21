#!/usr/bin/env python3
"""
STEP 8B read-only audit:
- compare final QM XYZ atom order with existing MOL/SDF topology carriers;
- verify formula, formal charge, stereochemical identity, and element sequence;
- identify which existing topology file can safely receive the final QM coordinates.

The script does not modify the repository.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import math
import sys

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

ROOT = Path("/data/BIOFILM_INSILICO_SUPPORT")

COMPOUNDS = {
    "compound_08": ROOT / "01_qm_orca_xtb/compound_08/geometries/compound_08_wb97xd3_def2svp_cpcm_water_opt_min.xyz",
    "compound_09": ROOT / "01_qm_orca_xtb/compound_09/geometries/compound_09_wb97xd3_def2svp_cpcm_water_minimum.xyz",
    "compound_11": ROOT / "01_qm_orca_xtb/compound_11/geometries/compound_11_wb97xd3_def2svp_cpcm_water_minimum.xyz",
}

SMILES_FILES = {
    name: ROOT / f"00_validated_inputs/ligands/{name}/{name}.smiles"
    for name in COMPOUNDS
}

THIOSEMICARBAZONE = Chem.MolFromSmarts("[C;X3](=[N;X2]-[N;X3]-[C;X3](=[S;X1])-[N;X3])")


def read_xyz(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 3:
        raise ValueError(f"Malformed XYZ: {path}")
    declared = int(lines[0].strip())
    atoms = []
    for line in lines[2:]:
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Malformed XYZ atom line in {path}: {line}")
        atoms.append((parts[0], tuple(float(x) for x in parts[1:4])))
    if declared != len(atoms):
        raise ValueError(f"{path}: declared {declared}, parsed {len(atoms)} atoms")
    return atoms


def load_single_molecule(path: Path):
    suffix = path.suffix.lower()
    if suffix == ".sdf":
        supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
        mols = [m for m in supplier if m is not None]
        if len(mols) != 1:
            return None, f"SDF molecules={len(mols)}"
        return mols[0], "OK"
    if suffix == ".mol":
        mol = Chem.MolFromMolFile(str(path), removeHs=False, sanitize=True)
        return (mol, "OK") if mol is not None else (None, "RDKit read failed")
    return None, "unsupported"


def no_h_canonical_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)), isomericSmiles=True, canonical=True)


def element_sequence(mol: Chem.Mol):
    return [atom.GetSymbol() for atom in mol.GetAtoms()]


def composition(seq):
    c = Counter(seq)
    order = ["C", "H", "N", "O", "S", "P", "F", "Cl", "Br", "I"]
    return " ".join(f"{e}={c[e]}" for e in order if c[e]) + "".join(
        f" {e}={c[e]}" for e in sorted(c) if e not in order
    )


def chiral_summary(mol: Chem.Mol):
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    return centers


def sequence_mismatch_positions(a, b, limit=12):
    out = []
    for i, (x, y) in enumerate(zip(a, b), start=1):
        if x != y:
            out.append((i, x, y))
            if len(out) >= limit:
                break
    if len(a) != len(b):
        out.append(("length", len(a), len(b)))
    return out


def print_thiosemicarbazone(mol: Chem.Mol):
    if THIOSEMICARBAZONE is None:
        print("  motif audit: SMARTS construction failed")
        return
    matches = mol.GetSubstructMatches(THIOSEMICARBAZONE)
    print(f"  thiosemicarbazone matches: {len(matches)}")
    for match in matches[:3]:
        atom_text = ", ".join(
            f"{idx + 1}:{mol.GetAtomWithIdx(idx).GetSymbol()}" for idx in match
        )
        bond_text = []
        for i, j in zip(match[:-1], match[1:]):
            bond = mol.GetBondBetweenAtoms(i, j)
            if bond is not None:
                bond_text.append(f"{i+1}-{j+1}:{bond.GetBondType()}")
        print(f"    atoms [{atom_text}]")
        print(f"    sequential bonds [{', '.join(bond_text)}]")


def main():
    print("=== STEP 8B READ-ONLY TOPOLOGY / ATOM-ORDER AUDIT ===")
    print(f"Repository: {ROOT}")
    print(f"RDKit: {Chem.rdBase.rdkitVersion}")

    for name, xyz_path in COMPOUNDS.items():
        print("\n" + "=" * 88)
        print(name)
        print("=" * 88)

        if not xyz_path.is_file():
            print(f"ERROR: missing final XYZ: {xyz_path}")
            continue

        xyz_atoms = read_xyz(xyz_path)
        xyz_seq = [e for e, _ in xyz_atoms]
        print(f"Final XYZ: {xyz_path.relative_to(ROOT)}")
        print(f"  atoms: {len(xyz_seq)}")
        print(f"  composition: {composition(xyz_seq)}")
        print(f"  element sequence head: {' '.join(xyz_seq[:20])}")
        print(f"  element sequence tail: {' '.join(xyz_seq[-20:])}")

        smiles = SMILES_FILES[name].read_text(encoding="utf-8").strip()
        validated = Chem.MolFromSmiles(smiles)
        if validated is None:
            print("ERROR: validated SMILES failed RDKit parsing")
            continue
        validated_h = Chem.AddHs(validated)
        validated_seq = element_sequence(validated_h)
        validated_canon = no_h_canonical_smiles(validated_h)
        print(f"Validated SMILES formula: {rdMolDescriptors.CalcMolFormula(validated_h)}")
        print(f"Validated formal charge: {Chem.GetFormalCharge(validated_h)}")
        print(f"Validated explicit-H atoms: {validated_h.GetNumAtoms()}")
        print(f"Validated chiral centers: {chiral_summary(validated_h)}")
        print(f"Validated SMILES-derived element sequence exact vs XYZ: {validated_seq == xyz_seq}")
        if validated_seq != xyz_seq:
            print(f"  first mismatches: {sequence_mismatch_positions(validated_seq, xyz_seq)}")
        if name in ("compound_08", "compound_09"):
            print_thiosemicarbazone(validated_h)

        geom_dir = xyz_path.parent
        candidates = sorted(
            p for p in geom_dir.rglob("*")
            if p.is_file() and p.suffix.lower() in {".sdf", ".mol"}
        )
        print(f"\nExisting MOL/SDF candidates found: {len(candidates)}")
        exact_candidates = []

        for path in candidates:
            mol, status = load_single_molecule(path)
            rel = path.relative_to(ROOT)
            if mol is None:
                print(f"  FAIL  {rel} | {status}")
                continue

            seq = element_sequence(mol)
            seq_exact = seq == xyz_seq
            formula = rdMolDescriptors.CalcMolFormula(mol)
            charge = Chem.GetFormalCharge(mol)
            canon = no_h_canonical_smiles(mol)
            identity_match = canon == validated_canon
            conf3d = bool(mol.GetNumConformers() and mol.GetConformer().Is3D())
            centers = chiral_summary(mol)

            mark = "CANDIDATE" if seq_exact and identity_match else "AUDIT"
            print(
                f"  {mark:9s} {rel}\n"
                f"    atoms={mol.GetNumAtoms()} bonds={mol.GetNumBonds()} "
                f"formula={formula} charge={charge} 3D={conf3d}\n"
                f"    element_sequence_exact={seq_exact} "
                f"stereochemical_identity_match={identity_match} "
                f"chiral_centers={centers}"
            )
            if not seq_exact:
                print(f"    first sequence mismatches: {sequence_mismatch_positions(seq, xyz_seq)}")
            if name in ("compound_08", "compound_09") and seq_exact and identity_match:
                print_thiosemicarbazone(mol)
            if seq_exact and identity_match:
                exact_candidates.append(path)

        print("\nRESULT")
        if exact_candidates:
            print("  Safe topology-carrier candidates (exact atom order + validated identity):")
            for p in exact_candidates:
                print(f"    - {p.relative_to(ROOT)}")
            print(
                "  These files are candidates for coordinate replacement with the final QM XYZ "
                "without distance-based bond perception."
            )
        else:
            print(
                "  WARNING: no existing MOL/SDF file simultaneously matched the final XYZ atom "
                "order and validated stereochemical identity."
            )

    print("\n=== AUDIT COMPLETE: repository was not modified ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
