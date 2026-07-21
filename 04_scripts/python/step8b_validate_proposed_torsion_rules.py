#!/usr/bin/env python3
"""
STEP 8B read-only validation of proposed ligand torsion rules.

Proposed chemistry-aware rules:
- compounds 08/09: rigidify the conjugated hydrazone N–N bond in
  C=N–NH–C(=S)–NH2;
- compound 11: rigidify the conjugated exocyclic NH2–thiadiazole bond;
- retain the non-conjugated thiadiazole-N–scaffold sigma bond as rotatable.

The script compares Meeko defaults with the proposed rules, maps BRANCH records
back to original MOL atom indices, and confirms that only the intended bond is removed.

The repository is not modified.
"""

from __future__ import annotations

from pathlib import Path
import re
import subprocess
import sys

from rdkit import Chem

ROOT = Path("/data/BIOFILM_INSILICO_SUPPORT")

CASES = {
    "compound_08": {
        "path": ROOT / "01_qm_orca_xtb/compound_08/geometries/compound_08_selected_for_xtb.mol",
        "smarts": "[C;X3](=[N;X2]-[N;X3]-[C;X3](=[S;X1])-[N;X3])",
        "bond_indices": ("2", "3"),
        "target_label": "conjugated hydrazone N-N",
        "target_smarts": Chem.MolFromSmarts("[C;X3](=[N;X2]-[N;X3]-[C;X3](=[S;X1])-[N;X3])"),
        "target_atoms_in_match": (1, 2),
    },
    "compound_09": {
        "path": ROOT / "01_qm_orca_xtb/compound_09/geometries/compound_09_rdkit_ETKDGv3_best.mol",
        "smarts": "[C;X3](=[N;X2]-[N;X3]-[C;X3](=[S;X1])-[N;X3])",
        "bond_indices": ("2", "3"),
        "target_label": "conjugated hydrazone N-N",
        "target_smarts": Chem.MolFromSmarts("[C;X3](=[N;X2]-[N;X3]-[C;X3](=[S;X1])-[N;X3])"),
        "target_atoms_in_match": (1, 2),
    },
    "compound_11": {
        "path": ROOT / "01_qm_orca_xtb/compound_11/geometries/compound_11_rdkit_ETKDGv3_best.mol",
        "smarts": "[N;X3;H2]-[c;R]",
        "bond_indices": ("1", "2"),
        "target_label": "conjugated exocyclic NH2-thiadiazole bond",
        "target_smarts": Chem.MolFromSmarts("[N;X3;H2]-[c;R]"),
        "target_atoms_in_match": (0, 1),
    },
}


def run_meeko(path: Path, rigid_smarts: str | None = None, bond_indices: tuple[str, str] | None = None) -> str:
    cmd = [
        "mk_prepare_ligand.py",
        "-i", str(path),
        "--charge_model", "gasteiger",
        "--add_index_map",
    ]
    if rigid_smarts is not None:
        cmd += ["-r", rigid_smarts, "-b", *bond_indices]
    cmd += ["-"]
    proc = subprocess.run(cmd, text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Meeko failed\nCOMMAND: {' '.join(cmd)}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    return proc.stdout


def parse_index_map(text: str) -> dict[int, int]:
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


def parse_branch_keys(text: str) -> set[tuple[int, int]]:
    mapping = parse_index_map(text)
    keys: set[tuple[int, int]] = set()
    for line in text.splitlines():
        fields = line.split()
        if fields and fields[0] == "BRANCH":
            a = mapping[int(fields[1])]
            b = mapping[int(fields[2])]
            keys.add(tuple(sorted((a, b))))
    return keys


def parse_torsdof(text: str) -> int:
    for line in text.splitlines():
        fields = line.split()
        if fields and fields[0] == "TORSDOF":
            return int(fields[1])
    raise RuntimeError("TORSDOF not found")


def target_bond_key(mol: Chem.Mol, query: Chem.Mol, atom_pair: tuple[int, int]) -> tuple[int, int]:
    matches = mol.GetSubstructMatches(query)
    if len(matches) != 1:
        raise RuntimeError(f"Expected one target SMARTS match, found {len(matches)}")
    match = matches[0]
    a = match[atom_pair[0]] + 1
    b = match[atom_pair[1]] + 1
    return tuple(sorted((a, b)))


def atom_label(mol: Chem.Mol, idx1: int) -> str:
    atom = mol.GetAtomWithIdx(idx1 - 1)
    return f"{idx1}:{atom.GetSymbol()}"


def describe_bond(mol: Chem.Mol, key: tuple[int, int]) -> str:
    a, b = key
    bond = mol.GetBondBetweenAtoms(a - 1, b - 1)
    if bond is None:
        return f"{a}-{b} NO_BOND"
    return (
        f"{atom_label(mol, a)}-{atom_label(mol, b)} "
        f"type={bond.GetBondType()} conjugated={bond.GetIsConjugated()} ring={bond.IsInRing()}"
    )


def main() -> int:
    print("=== STEP 8B: VALIDATION OF PROPOSED TORSION RULES ===")
    print(f"RDKit: {Chem.rdBase.rdkitVersion}")
    print("Mode: read-only")

    all_ok = True

    for name, cfg in CASES.items():
        print("\n" + "=" * 96)
        print(name)
        print("=" * 96)

        path = cfg["path"]
        mol = Chem.MolFromMolFile(str(path), removeHs=False, sanitize=True)
        if mol is None:
            raise RuntimeError(f"RDKit failed to read {path}")

        default_pdbqt = run_meeko(path)
        proposed_pdbqt = run_meeko(path, cfg["smarts"], cfg["bond_indices"])

        default_keys = parse_branch_keys(default_pdbqt)
        proposed_keys = parse_branch_keys(proposed_pdbqt)
        removed = default_keys - proposed_keys
        added = proposed_keys - default_keys
        target = target_bond_key(
            mol,
            cfg["target_smarts"],
            cfg["target_atoms_in_match"],
        )

        default_tors = parse_torsdof(default_pdbqt)
        proposed_tors = parse_torsdof(proposed_pdbqt)

        print(f"Topology file: {path.relative_to(ROOT)}")
        print(f"Rigidify SMARTS: {cfg['smarts']}")
        print(f"SMARTS bond indices: {' '.join(cfg['bond_indices'])}")
        print(f"Target: {cfg['target_label']}")
        print(f"Target bond: {describe_bond(mol, target)}")
        print(f"Default TORSDOF:  {default_tors}")
        print(f"Proposed TORSDOF: {proposed_tors}")
        print(f"Removed branches: {len(removed)}")
        for key in sorted(removed):
            print(f"  REMOVED: {describe_bond(mol, key)}")
        print(f"Added branches: {len(added)}")
        for key in sorted(added):
            print(f"  ADDED:   {describe_bond(mol, key)}")

        target_default = target in default_keys
        target_proposed = target in proposed_keys
        only_target_removed = removed == {target}
        no_added = not added
        torsion_decrement_one = proposed_tors == default_tors - 1

        print(f"Target was rotatable by default: {target_default}")
        print(f"Target remains rotatable after rule: {target_proposed}")
        print(f"Only intended target removed: {only_target_removed}")
        print(f"No new branches introduced: {no_added}")
        print(f"TORSDOF decreased by exactly 1: {torsion_decrement_one}")

        case_ok = (
            target_default
            and not target_proposed
            and only_target_removed
            and no_added
            and torsion_decrement_one
        )
        print(f"CASE STATUS: {'PASS' if case_ok else 'FAIL'}")
        all_ok = all_ok and case_ok

        if name == "compound_11":
            ring_n_scaffold = Chem.MolFromSmarts("[n;R]-[C;X4]")
            matches = mol.GetSubstructMatches(ring_n_scaffold)
            if len(matches) == 1:
                a, b = matches[0]
                key = tuple(sorted((a + 1, b + 1)))
                bond = mol.GetBondBetweenAtoms(a, b)
                print(
                    "Retained genuine sigma rotor (ring N-scaffold): "
                    f"{describe_bond(mol, key)} "
                    f"rotatable_after_rule={key in proposed_keys}"
                )
            else:
                print(f"WARNING: ring N-scaffold matches={len(matches)}")

    print("\n" + "=" * 96)
    print(f"OVERALL STATUS: {'PASS' if all_ok else 'FAIL'}")
    print("Repository was not modified.")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
