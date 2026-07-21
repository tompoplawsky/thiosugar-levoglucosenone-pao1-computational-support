from pathlib import Path
from tempfile import TemporaryDirectory
from subprocess import run, PIPE
from hashlib import sha256
from datetime import datetime, timezone
from collections import Counter
from math import dist

import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from scipy.optimize import linear_sum_assignment

root = Path("/data/BIOFILM_INSILICO_SUPPORT")

receptor = root / (
    "00_validated_inputs/receptors/LasR_2UV0/"
    "LasR_2UV0_prepared_chainF.pdb"
)

source_pdbqt = root / (
    "03_docking/vina_outputs/production/LasR_2UV0_2026-07-16/"
    "selected_representative/"
    "compound_11_LasR_2UV0_consensus_cluster1_run01_mode06.pdbqt"
)

ligand_sdf = root / (
    "03_docking/plip_inputs/LasR_2UV0_2026-07-17/"
    "meeko_export_audit/"
    "compound_11_LasR_2UV0_consensus_run01_mode06_meeko.sdf"
)

outdir = root / "03_docking/plip_inputs/LasR_2UV0_2026-07-17/final"
complex_out = outdir / "LasR_2UV0_compound_11_consensus_complex_FINAL.pdb"
record_out = outdir / "LasR_2UV0_compound_11_consensus_complex_FINAL_record.txt"

for path in (complex_out, record_out):
    if path.exists():
        raise SystemExit(f"REFUSING TO OVERWRITE: {path}")

for path in (receptor, source_pdbqt, ligand_sdf):
    if not path.is_file():
        raise SystemExit(f"MISSING INPUT: {path}")

outdir.mkdir(parents=True, exist_ok=True)


def digest(path: Path) -> str:
    h = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def xyz(line: str):
    return np.array(
        [
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        ],
        dtype=float,
    )


def parse_conect(line: str):
    return [
        int(line[i:i + 5])
        for i in range(6, len(line), 5)
        if line[i:i + 5].strip()
    ]


supplier = Chem.SDMolSupplier(
    str(ligand_sdf),
    removeHs=False,
    sanitize=True,
)
molecules = [mol for mol in supplier if mol is not None]

if len(molecules) != 1:
    raise SystemExit(
        f"Expected one valid molecule in SDF; found {len(molecules)}"
    )

mol = molecules[0]
formula = rdMolDescriptors.CalcMolFormula(mol)
formal_charge = Chem.GetFormalCharge(mol)
heavy_atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
elements = Counter(a.GetSymbol() for a in heavy_atoms)

if formula != "C8H9N3O3S2":
    raise SystemExit(f"Unexpected formula: {formula}")
if formal_charge != 0:
    raise SystemExit(f"Unexpected formal charge: {formal_charge}")
if len(heavy_atoms) != 16:
    raise SystemExit(f"Unexpected heavy-atom count: {len(heavy_atoms)}")
if elements != Counter({"C": 8, "N": 3, "O": 3, "S": 2}):
    raise SystemExit(f"Unexpected element counts: {dict(elements)}")

sdf_conf = mol.GetConformer()
sdf_atoms = [
    (
        atom.GetSymbol(),
        np.array(sdf_conf.GetAtomPosition(atom.GetIdx()), dtype=float),
    )
    for atom in mol.GetAtoms()
]

obabel_version = run(
    [
        "conda", "run", "-n", "plip_env",
        "obabel", "-V",
    ],
    check=True,
    text=True,
    stdout=PIPE,
    stderr=PIPE,
).stdout.strip()

with TemporaryDirectory() as tmpdir:
    ligand_tmp = Path(tmpdir) / "compound_11_meeko_export.pdb"

    conversion = run(
        [
            "conda", "run", "-n", "plip_env",
            "obabel",
            str(ligand_sdf),
            "-O", str(ligand_tmp),
        ],
        check=True,
        text=True,
        stdout=PIPE,
        stderr=PIPE,
    )

    tmp_lines = ligand_tmp.read_text().splitlines()
    ligand_atoms_old = [
        line for line in tmp_lines
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    conect_old = [
        line for line in tmp_lines
        if line.startswith("CONECT")
    ]

if len(ligand_atoms_old) != mol.GetNumAtoms():
    raise SystemExit(
        f"Open Babel atom count mismatch: "
        f"{len(ligand_atoms_old)} versus {mol.GetNumAtoms()}"
    )

if not conect_old:
    raise SystemExit("Open Babel generated no ligand CONECT records")

pdb_atoms = []
for line in ligand_atoms_old:
    element = line[76:78].strip()
    if not element:
        atom_name = line[12:16].strip()
        element = "".join(c for c in atom_name if c.isalpha())[:1].upper()

    pdb_atoms.append((element, xyz(line)))

if Counter(el for el, _ in pdb_atoms) != Counter(el for el, _ in sdf_atoms):
    raise SystemExit("Element counts changed during SDF-to-PDB conversion")

matched_distances = []
for element in sorted(set(el for el, _ in sdf_atoms)):
    a = np.array([coord for el, coord in sdf_atoms if el == element])
    b = np.array([coord for el, coord in pdb_atoms if el == element])
    cost = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    rows, cols = linear_sum_assignment(cost)
    matched_distances.extend(cost[rows, cols].tolist())

max_conversion_shift = max(matched_distances)
rms_conversion_shift = float(
    np.sqrt(np.mean(np.square(matched_distances)))
)

if max_conversion_shift > 0.002:
    raise SystemExit(
        f"SDF-to-PDB coordinate shift too large: "
        f"{max_conversion_shift:.6f} A"
    )

receptor_atoms = [
    line.rstrip()
    for line in receptor.read_text().splitlines()
    if line.startswith(("ATOM  ", "HETATM"))
]

if len(receptor_atoms) != 1295:
    raise SystemExit(
        f"Unexpected receptor atom count: {len(receptor_atoms)}"
    )

max_receptor_serial = max(int(line[6:11]) for line in receptor_atoms)
receptor_ter_serial = max_receptor_serial + 1
first_ligand_serial = receptor_ter_serial + 1

old_serials = [int(line[6:11]) for line in ligand_atoms_old]
serial_map = {
    old: first_ligand_serial + index
    for index, old in enumerate(old_serials)
}

ligand_records = []

for index, line in enumerate(ligand_atoms_old):
    serial = first_ligand_serial + index
    atom_name = line[12:16].strip()
    element = line[76:78].strip()

    if not element:
        element = "".join(
            c for c in atom_name if c.isalpha()
        )[:1].upper()

    x, y, z = xyz(line)

    ligand_records.append(
        f"HETATM{serial:5d} {atom_name:>4s} "
        f"L11 Z{901:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.00:6.2f}{0.00:6.2f}"
        f"          {element:>2s}"
    )

mapped_conect = []

for line in conect_old:
    values = parse_conect(line)

    if not values:
        continue

    if any(value not in serial_map for value in values):
        raise SystemExit(f"Unknown atom serial in CONECT: {line}")

    mapped_conect.append(
        "CONECT" +
        "".join(f"{serial_map[value]:5d}" for value in values)
    )

last_ligand_serial = first_ligand_serial + len(ligand_records) - 1
ligand_ter_serial = last_ligand_serial + 1

complex_lines = [
    "HEADER    LASR-COMPOUND 11 VALIDATED DOCKING COMPLEX",
    "REMARK 900 RECEPTOR LASR 2UV0 PREPARED CHAIN F",
    "REMARK 900 LIGAND COMPOUND 11 RESNAME L11 CHAIN Z RESIDUE 901",
    "REMARK 900 REPRESENTATIVE RUN01 SEED 2026071661 MODE 06",
    "REMARK 900 VINARDO SCORE -3.5860 KCAL/MOL",
    "REMARK 900 LIGAND CHEMISTRY RECONSTRUCTED BY MEEKO 0.7.1",
    "REMARK 900 LIGAND CONNECTIVITY TRANSFERRED FROM VALIDATED SDF",
    *receptor_atoms,
    f"TER   {receptor_ter_serial:5d}      GLU F 168",
    *ligand_records,
    f"TER   {ligand_ter_serial:5d}      L11 Z 901",
    *mapped_conect,
    "END",
]

complex_out.write_text("\n".join(complex_lines) + "\n")

coordinate_serials = [
    int(line[6:11])
    for line in complex_lines
    if line.startswith(("ATOM  ", "HETATM", "TER   "))
]

serial_counts = Counter(coordinate_serials)
duplicates = sorted(
    serial for serial, count in serial_counts.items()
    if count > 1
)

if duplicates:
    raise SystemExit(f"Duplicate serial numbers: {duplicates}")

final_ligand_xyz = [
    xyz(line)
    for line in ligand_records
]

minimum_contact = min(
    dist(
        xyz(receptor_line).tolist(),
        ligand_coord.tolist(),
    )
    for receptor_line in receptor_atoms
    for ligand_coord in final_ligand_xyz
    if receptor_line[76:78].strip().upper() != "H"
)

record_text = f"""LasR–compound_11 final PLIP complex construction record

Date UTC: {datetime.now(timezone.utc).isoformat()}

Inputs:
- receptor: {receptor}
- selected docking pose: {source_pdbqt}
- Meeko-reconstructed ligand: {ligand_sdf}

Output:
- final PLIP complex: {complex_out}

Ligand validation:
- formula: {formula}
- formal charge: {formal_charge}
- total atoms: {mol.GetNumAtoms()}
- heavy atoms: {len(heavy_atoms)}
- heavy-atom elements: {dict(sorted(elements.items()))}

Coordinate validation:
- maximum SDF-to-PDB atom displacement: {max_conversion_shift:.6f} A
- RMS SDF-to-PDB atom displacement: {rms_conversion_shift:.6f} A
- minimum receptor-ligand atom distance: {minimum_contact:.4f} A

PDB validation:
- receptor atoms: {len(receptor_atoms)}
- receptor TER serial: {receptor_ter_serial}
- ligand atoms: {len(ligand_records)}
- ligand atom serial range: {first_ligand_serial}-{last_ligand_serial}
- ligand TER serial: {ligand_ter_serial}
- ligand CONECT records: {len(mapped_conect)}
- duplicate coordinate/TER serials: none

Software:
- RDKit: {rdkit.__version__}
- Open Babel: {obabel_version}
- Meeko export source: Meeko 0.7.1

SHA-256:
- receptor: {digest(receptor)}
- selected PDBQT: {digest(source_pdbqt)}
- ligand SDF: {digest(ligand_sdf)}
- final complex: {digest(complex_out)}

Status: PASS
"""

record_out.write_text(record_text)

print(record_text)
print("===== FINAL COMPLEX LAST 15 LINES =====")
for line in complex_lines[-15:]:
    print(line)

print(f"\nFINAL COMPLEX: {complex_out}")
print(f"RECORD:        {record_out}")
