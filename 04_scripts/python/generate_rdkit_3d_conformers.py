#!/usr/bin/env python3

from pathlib import Path
import argparse
import csv
import json
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


EXPECTED = {
    "compound_08": {
        "smiles_file": "00_validated_inputs/ligands/compound_08/compound_08.smiles",
        "out_dir": "01_qm_orca_xtb/compound_08/geometries",
        "formula": "C33H45N3O19S2",
        "mw": 851.85,
        "nconf": 250,
    },
    "compound_09": {
        "smiles_file": "00_validated_inputs/ligands/compound_09/compound_09.smiles",
        "out_dir": "01_qm_orca_xtb/compound_09/geometries",
        "formula": "C18H25N3O9S2",
        "mw": 491.53,
        "nconf": 150,
    },
    "compound_11": {
        "smiles_file": "00_validated_inputs/ligands/compound_11/compound_11.smiles",
        "out_dir": "01_qm_orca_xtb/compound_11/geometries",
        "formula": "C8H9N3O3S2",
        "mw": 259.30,
        "nconf": 75,
    },
}


THIOSEMICARBAZONE_SMARTS = Chem.MolFromSmarts("[#6]=[#7]-[#7]-[#6](=[#16])-[#7]")


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def read_smiles(path: Path) -> str:
    if not path.exists():
        fail(f"SMILES file not found: {path}")
    txt = path.read_text().strip()
    if not txt:
        fail(f"Empty SMILES file: {path}")
    return txt.split()[0]


def check_basic_identity(mol: Chem.Mol, compound: str, expected: dict) -> dict:
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mw = Descriptors.MolWt(mol)
    formal_charge = Chem.GetFormalCharge(mol)
    chiral = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)

    if formula != expected["formula"]:
        fail(f"{compound}: formula mismatch: parsed {formula}, expected {expected['formula']}")

    if abs(mw - expected["mw"]) > 0.15:
        fail(f"{compound}: MW mismatch: parsed {mw:.2f}, expected {expected['mw']:.2f}")

    if formal_charge != 0:
        fail(f"{compound}: non-zero formal charge detected: {formal_charge}")

    if compound in ("compound_08", "compound_09"):
        if not mol.HasSubstructMatch(THIOSEMICARBAZONE_SMARTS):
            fail(f"{compound}: thiosemicarbazone motif C=N-N-C(=S)-N not detected")

    return {
        "formula": formula,
        "mw": round(mw, 4),
        "formal_charge": formal_charge,
        "chiral_centers": chiral,
        "num_atoms_noH": mol.GetNumAtoms(),
    }


def embed_conformers(molH: Chem.Mol, nconf: int, seed: int):
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.enforceChirality = True
    params.useSmallRingTorsions = True
    params.useBasicKnowledge = True
    params.useExpTorsionAnglePrefs = True
    params.pruneRmsThresh = 0.5
    params.numThreads = 0

    conf_ids = list(AllChem.EmbedMultipleConfs(molH, numConfs=nconf, params=params))

    if len(conf_ids) == 0:
        params.useRandomCoords = True
        params.boxSizeMult = 4.0
        conf_ids = list(AllChem.EmbedMultipleConfs(molH, numConfs=nconf, params=params))

    if len(conf_ids) == 0:
        fail("RDKit embedding failed: no conformers generated")

    return conf_ids


def optimize_conformers(molH: Chem.Mol, conf_ids):
    if AllChem.MMFFHasAllMoleculeParams(molH):
        ff_name = "MMFF94s"
        results = AllChem.MMFFOptimizeMoleculeConfs(
            molH,
            numThreads=0,
            maxIters=2000,
            mmffVariant="MMFF94s",
        )
    else:
        ff_name = "UFF"
        results = AllChem.UFFOptimizeMoleculeConfs(
            molH,
            numThreads=0,
            maxIters=2000,
        )

    if len(results) != len(conf_ids):
        fail("Optimization result count does not match conformer count")

    ranked = []
    for cid, result in zip(conf_ids, results):
        not_converged, energy = result
        ranked.append({
            "conf_id": int(cid),
            "not_converged": int(not_converged),
            "energy": float(energy),
        })

    ranked.sort(key=lambda x: x["energy"])
    return ff_name, ranked


def write_outputs(molH: Chem.Mol, compound: str, smiles: str, out_dir: Path, ff_name: str, ranked: list, identity: dict):
    out_dir.mkdir(parents=True, exist_ok=True)

    all_sdf = out_dir / f"{compound}_rdkit_ETKDGv3_all_conformers.sdf"
    best_sdf = out_dir / f"{compound}_rdkit_ETKDGv3_best.sdf"
    best_mol = out_dir / f"{compound}_rdkit_ETKDGv3_best.mol"
    best_xyz = out_dir / f"{compound}_rdkit_ETKDGv3_best.xyz"
    report_csv = out_dir / f"{compound}_rdkit_ETKDGv3_conformer_report.csv"
    report_json = out_dir / f"{compound}_rdkit_ETKDGv3_generation_report.json"

    writer = Chem.SDWriter(str(all_sdf))
    for rank, row in enumerate(ranked, start=1):
        cid = row["conf_id"]
        molH.SetProp("_Name", f"{compound}_conf_{rank:03d}_cid_{cid}")
        molH.SetProp("compound", compound)
        molH.SetProp("source_smiles", smiles)
        molH.SetProp("rank", str(rank))
        molH.SetProp("conf_id", str(cid))
        molH.SetProp("force_field", ff_name)
        molH.SetProp("energy", f"{row['energy']:.8f}")
        molH.SetProp("not_converged", str(row["not_converged"]))
        molH.SetProp("rdkit_embedding", "ETKDGv3")
        writer.write(molH, confId=cid)
    writer.close()

    best_cid = ranked[0]["conf_id"]
    writer = Chem.SDWriter(str(best_sdf))
    molH.SetProp("_Name", f"{compound}_best_cid_{best_cid}")
    molH.SetProp("compound", compound)
    molH.SetProp("source_smiles", smiles)
    molH.SetProp("rank", "1")
    molH.SetProp("conf_id", str(best_cid))
    molH.SetProp("force_field", ff_name)
    molH.SetProp("energy", f"{ranked[0]['energy']:.8f}")
    molH.SetProp("not_converged", str(ranked[0]["not_converged"]))
    molH.SetProp("rdkit_embedding", "ETKDGv3")
    writer.write(molH, confId=best_cid)
    writer.close()

    Chem.MolToMolFile(molH, str(best_mol), confId=best_cid)
    Chem.MolToXYZFile(molH, str(best_xyz), confId=best_cid)

    with report_csv.open("w", newline="") as f:
        fieldnames = ["rank", "conf_id", "force_field", "energy", "not_converged"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for rank, row in enumerate(ranked, start=1):
            w.writerow({
                "rank": rank,
                "conf_id": row["conf_id"],
                "force_field": ff_name,
                "energy": f"{row['energy']:.8f}",
                "not_converged": row["not_converged"],
            })

    metadata = {
        "compound": compound,
        "source_smiles": smiles,
        "rdkit_embedding": "ETKDGv3",
        "force_field": ff_name,
        "n_conformers_generated": len(ranked),
        "best_conf_id": best_cid,
        "best_energy": ranked[0]["energy"],
        "identity_check": identity,
        "outputs": {
            "all_conformers_sdf": str(all_sdf),
            "best_sdf": str(best_sdf),
            "best_mol": str(best_mol),
            "best_xyz": str(best_xyz),
            "report_csv": str(report_csv),
        },
        "warning": "MMFF/UFF geometry is only a pre-QM starting structure, not a final quantum-chemical geometry.",
    }

    report_json.write_text(json.dumps(metadata, indent=2))

    print(f"OK: {compound}")
    print(f"Generated conformers: {len(ranked)}")
    print(f"Force field: {ff_name}")
    print(f"Best conformer ID: {best_cid}")
    print(f"Best energy: {ranked[0]['energy']:.8f}")
    print(f"Best XYZ: {best_xyz}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".", help="BIOFILM_INSILICO_SUPPORT root directory")
    parser.add_argument("--compound", required=True, choices=EXPECTED.keys())
    parser.add_argument("--seed", type=int, default=20260710)
    args = parser.parse_args()

    root = Path(args.root).resolve()
    cfg = EXPECTED[args.compound]

    smiles_path = root / cfg["smiles_file"]
    out_dir = root / cfg["out_dir"]

    smiles = read_smiles(smiles_path)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        fail(f"{args.compound}: RDKit could not parse SMILES")

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    identity = check_basic_identity(mol, args.compound, cfg)

    molH = Chem.AddHs(mol)
    conf_ids = embed_conformers(molH, cfg["nconf"], args.seed)
    ff_name, ranked = optimize_conformers(molH, conf_ids)

    write_outputs(molH, args.compound, smiles, out_dir, ff_name, ranked, identity)


if __name__ == "__main__":
    main()
