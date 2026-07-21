#!/usr/bin/env python3
"""
STEP 8B read-only audit of Meeko default rotatable bonds.

For each ligand:
- runs mk_prepare_ligand.py to stdout;
- parses REMARK INDEX MAP and BRANCH records;
- maps every PDBQT branch back to the original MOL atom indices;
- reports bond type, atom environments, ring status, and conjugation;
- explicitly audits the thiosemicarbazone bonds in compounds 08/09;
- explicitly audits the exocyclic amino and ring-N–scaffold bonds in compound 11.

The repository is not modified.
"""

from __future__ import annotations

from pathlib import Path
import re
import subprocess
import sys

from rdkit import Chem

ROOT = Path("/data/BIOFILM_INSILICO_SUPPORT")

LIGANDS = {
    "compound_08": ROOT / "01_qm_orca_xtb/compound_08/geometries/compound_08_selected_for_xtb.mol",
    "compound_09": ROOT / "01_qm_orca_xtb/compound_09/geometries/compound_09_rdkit_ETKDGv3_best.mol",
    "compound_11": ROOT / "01_qm_orca_xtb/compound_11/geometries/compound_11_rdkit_ETKDGv3_best.mol",
}

THIOSEMICARBAZONE = Chem.MolFromSmarts("[C;X3](=[N;X2]-[N;X3]-[C;X3](=[S;X1])-[N;X3])")
EXOCYCLIC_AMINO_11 = Chem.MolFromSmarts("[N;X3;H2]-[c;R]")
RING_N_SCAFFOLD_11 = Chem.MolFromSmarts("[n;R]-[C;X4]")


def atom_label(atom: Chem.Atom) -> str:
    idx = atom.GetIdx() + 1
    h = atom.GetTotalNumHs(includeNeighbors=True)
    charge = atom.GetFormalCharge()
    aromatic = "arom" if atom.GetIsAromatic() else "aliph"
    return f"{idx}:{atom.GetSymbol()}(H{h},q{charge},{aromatic})"


def neighbor_summary(atom: Chem.Atom) -> str:
    parts = []
    for nbr in atom.GetNeighbors():
        bond = atom.GetOwningMol().GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
        parts.append(f"{nbr.GetIdx()+1}:{nbr.GetSymbol()}[{bond.GetBondType()}]")
    return ", ".join(parts)


def parse_index_map(text: str) -> dict[int, int]:
    """Return PDBQT serial -> original MOL 1-based index."""
    mapping: dict[int, int] = {}
    for line in text.splitlines():
        if not line.startswith("REMARK INDEX MAP"):
            continue
        nums = [int(x) for x in re.findall(r"\d+", line[len("REMARK INDEX MAP"):])]
        if len(nums) % 2:
            raise RuntimeError(f"Odd INDEX MAP field count: {line}")
        for original_idx, pdbqt_serial in zip(nums[0::2], nums[1::2]):
            mapping[pdbqt_serial] = original_idx
    return mapping


def parse_branches(text: str) -> list[tuple[int, int]]:
    branches = []
    for line in text.splitlines():
        fields = line.split()
        if fields and fields[0] == "BRANCH":
            branches.append((int(fields[1]), int(fields[2])))
    return branches


def run_meeko(path: Path) -> str:
    cmd = [
        "mk_prepare_ligand.py",
        "-i", str(path),
        "--charge_model", "gasteiger",
        "--add_index_map",
        "-"
    ]
    proc = subprocess.run(cmd, text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Meeko failed for {path}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return proc.stdout


def bond_key(i: int, j: int) -> tuple[int, int]:
    return tuple(sorted((i, j)))


def report_special_bonds(name: str, mol: Chem.Mol, branch_keys: set[tuple[int, int]]) -> None:
    print("\nSPECIAL CONJUGATED-BOND AUDIT")

    if name in {"compound_08", "compound_09"}:
        matches = mol.GetSubstructMatches(THIOSEMICARBAZONE)
        print(f"  thiosemicarbazone matches: {len(matches)}")
        for match in matches:
            # SMARTS atom order: C(imine), N(imine), N(hydrazinic), C(thioamide), S, N(amino)
            ci, ni, nh, ct, s, na = match
            checks = [
                ("C=N imine", ci, ni),
                ("N–N hydrazone", ni, nh),
                ("N–C(=S) thioamide", nh, ct),
                ("C=S", ct, s),
                ("C–NH2 thioamide", ct, na),
            ]
            for label, a, b in checks:
                bond = mol.GetBondBetweenAtoms(a, b)
                key = bond_key(a + 1, b + 1)
                print(
                    f"  {label:23s} {a+1}-{b+1} "
                    f"{mol.GetAtomWithIdx(a).GetSymbol()}-{mol.GetAtomWithIdx(b).GetSymbol()} "
                    f"type={bond.GetBondType()} conjugated={bond.GetIsConjugated()} "
                    f"Meeko_branch={'YES' if key in branch_keys else 'NO'}"
                )

    if name == "compound_11":
        amino_matches = mol.GetSubstructMatches(EXOCYCLIC_AMINO_11)
        ring_n_matches = mol.GetSubstructMatches(RING_N_SCAFFOLD_11)
        print(f"  exocyclic amino-to-ring matches: {len(amino_matches)}")
        for n, c in amino_matches:
            bond = mol.GetBondBetweenAtoms(n, c)
            key = bond_key(n + 1, c + 1)
            print(
                f"  exocyclic NH2–ring C     {n+1}-{c+1} "
                f"type={bond.GetBondType()} conjugated={bond.GetIsConjugated()} "
                f"Meeko_branch={'YES' if key in branch_keys else 'NO'}"
            )
        print(f"  ring-N-to-scaffold matches: {len(ring_n_matches)}")
        for n, c in ring_n_matches:
            bond = mol.GetBondBetweenAtoms(n, c)
            key = bond_key(n + 1, c + 1)
            print(
                f"  ring N–scaffold C        {n+1}-{c+1} "
                f"type={bond.GetBondType()} conjugated={bond.GetIsConjugated()} "
                f"Meeko_branch={'YES' if key in branch_keys else 'NO'}"
            )


def main() -> int:
    print("=== STEP 8B: MEEKO DEFAULT ROTATABLE-BOND MAPPING ===")
    print(f"RDKit: {Chem.rdBase.rdkitVersion}")
    print("Mode: read-only")

    for name, path in LIGANDS.items():
        print("\n" + "=" * 100)
        print(name)
        print("=" * 100)
        mol = Chem.MolFromMolFile(str(path), removeHs=False, sanitize=True)
        if mol is None:
            raise RuntimeError(f"RDKit could not read {path}")

        pdbqt = run_meeko(path)
        index_map = parse_index_map(pdbqt)
        branches = parse_branches(pdbqt)

        print(f"Topology file: {path.relative_to(ROOT)}")
        print(f"Atoms={mol.GetNumAtoms()} Bonds={mol.GetNumBonds()} Meeko branches={len(branches)}")
        print(f"Mapped PDBQT atoms={len(index_map)}")

        branch_keys: set[tuple[int, int]] = set()
        print("\nDEFAULT MEEKO BRANCHES MAPPED TO ORIGINAL MOL")
        for number, (pdb_a, pdb_b) in enumerate(branches, start=1):
            if pdb_a not in index_map or pdb_b not in index_map:
                print(f"  {number:2d}. PDBQT {pdb_a}-{pdb_b}: mapping missing")
                continue
            orig_a = index_map[pdb_a]
            orig_b = index_map[pdb_b]
            branch_keys.add(bond_key(orig_a, orig_b))
            a = mol.GetAtomWithIdx(orig_a - 1)
            b = mol.GetAtomWithIdx(orig_b - 1)
            bond = mol.GetBondBetweenAtoms(orig_a - 1, orig_b - 1)
            if bond is None:
                print(f"  {number:2d}. PDBQT {pdb_a}-{pdb_b} -> MOL {orig_a}-{orig_b}: NO BOND")
                continue

            print(
                f"  {number:2d}. PDBQT {pdb_a:>2}-{pdb_b:<2} -> MOL {orig_a:>2}-{orig_b:<2} | "
                f"{atom_label(a)} -- {atom_label(b)} | "
                f"type={bond.GetBondType()} ring={bond.IsInRing()} "
                f"conjugated={bond.GetIsConjugated()} stereo={bond.GetStereo()}"
            )
            print(f"      atom {orig_a} neighbors: {neighbor_summary(a)}")
            print(f"      atom {orig_b} neighbors: {neighbor_summary(b)}")

        report_special_bonds(name, mol, branch_keys)

    print("\n=== AUDIT COMPLETE: repository was not modified ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
