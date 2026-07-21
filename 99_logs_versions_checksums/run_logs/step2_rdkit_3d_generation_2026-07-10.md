# Step 2 — RDKit 3D ligand geometry generation

Date: 2026-07-10
Workspace: /data/BIOFILM_INSILICO_SUPPORT/

Software:
- RDKit 2025.09.5
- Open Babel 3.1.0

Method:
- 3D conformers generated from validated SMILES using RDKit ETKDGv3.
- Explicit hydrogens and stereochemistry enforcement were used.
- Initial minimization used MMFF94s.
- MMFF94s geometries are pre-QM starting structures only.

Generated conformers:
- compound_08: 116 conformers retained; best conformer ID 104; best energy 71.03987918
- compound_09: 95 conformers retained; best conformer ID 14; best energy 51.24137133
- compound_11: 15 conformers retained; best conformer ID 11; best energy 16.47833698

Automated checks completed:
- RDKit reads all best SDF files.
- Open Babel converts all best SDF files.
- Formula and MW match expected validated ligand identities.
- Formal charge is zero for all ligands.
- Thiosemicarbazone motif C=N-N-C(=S)-N detected in compound_08 and compound_09.
- Canonical isomeric SMILES from best SDF matches canonical isomeric SMILES from validated input SMILES for all ligands.

Current status:
- Technical 3D generation completed.
- Step 2 is not yet closed.
- Pending: visual inspection of ring conformations, levoglucosenone core, hydrazone C=N geometry, N-N connectivity, and close contacts.
- Do not proceed to xTB or ORCA before visual inspection.

## Visual 3D inspection and conformer selection

Visual inspection was performed in PyMOL for compound_11, compound_09, and compound_08.

QC PNG files were generated for the inspected structures:
- compound_11: 01_qm_orca_xtb/compound_11/geometries/compound_11_rdkit_ETKDGv3_visual_qc.png
- compound_09: 01_qm_orca_xtb/compound_09/geometries/compound_09_rdkit_ETKDGv3_visual_qc.png
- compound_08: 01_qm_orca_xtb/compound_08/geometries/compound_08_rdkit_ETKDGv3_visual_qc.png
- compound_08 top10 panel: 01_qm_orca_xtb/compound_08/geometries/compound_08_top10_qc/compound_08_top10_visual_qc.png

compound_11 and compound_09 best MMFF94s conformers were visually accepted as chemically reasonable pre-QM starting structures.

For compound_08, the lowest-energy MMFF94s conformer was visually compact/folded. Because compound_08 is a large flexible peracetylated disaccharide derivative, the top 10 MMFF94s-ranked conformers were inspected. compound_08_top10_2.mol was selected as the xTB starting geometry because it appeared less artificially compact while preserving normal ring geometry and thiosemicarbazone connectivity. compound_08_top10_4.mol was retained as backup.

Selected compound_08 xTB input geometry:
- 01_qm_orca_xtb/compound_08/geometries/compound_08_selected_for_xtb.mol
- 01_qm_orca_xtb/compound_08/geometries/compound_08_selected_for_xtb.xyz

Backup compound_08 geometry:
- 01_qm_orca_xtb/compound_08/geometries/compound_08_backup_conformer_top4.mol
- 01_qm_orca_xtb/compound_08/geometries/compound_08_backup_conformer_top4.xyz

## Step 2 closure

Step 2 status: CLOSED.

Validated 3D ligand geometries were generated from final SMILES and inspected before quantum-chemical calculations.

Accepted pre-QM starting geometries:
- compound_11: 01_qm_orca_xtb/compound_11/geometries/compound_11_rdkit_ETKDGv3_best.xyz
- compound_09: 01_qm_orca_xtb/compound_09/geometries/compound_09_rdkit_ETKDGv3_best.xyz
- compound_08: 01_qm_orca_xtb/compound_08/geometries/compound_08_selected_for_xtb.xyz

Next pipeline step: Step 3 — xTB / ORCA quantum-chemical calculations.
