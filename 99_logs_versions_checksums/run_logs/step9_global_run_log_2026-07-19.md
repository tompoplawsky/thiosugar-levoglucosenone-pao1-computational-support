# STEP 9 global run log

**Repository:** `/data/BIOFILM_INSILICO_SUPPORT`  
**Closure generated UTC:** `2026-07-19T17:37:21.483450+00:00`

## Scope

STEP 9 covered receptor-specific redocking, LasR production docking, exploratory
PqsR docking, spatial classification, clustering/representative selection, PLIP
interaction profiling, supervised ChimeraX visual QC, checksums, and verified backups.

## Frozen docking protocol

- AutoDock Vina 1.2.5
- scoring function: Vinardo
- CPU: 8
- exhaustiveness: 64
- num_modes: 20
- energy_range: 5 kcal/mol
- receptor-specific 30 x 30 x 30 A boxes
- five recorded seeds for each production/exploratory receptor-ligand pair

## Redocking

- LasR/2UV0 reference ligand: OHN
- LasR validation: PASS
- best heavy-atom RMSD: 0.6911 A
- centroid displacement: 0.3053 A
- PqsR/4JVI reference ligand: QZN
- PqsR validation: FAIL after documented V1, V2, V3, and V4 workflows
- PqsR production poses are therefore exploratory structural rationalizations only

Redocking run-record files detected: 17

## LasR production docking

- compounds: 08, 09, 11
- independent runs: 5 per compound
- run-record files detected: 15
- 100 poses per compound

Spatial classification:
- compound 08: 0/100 orthosteric; 0/100 partial-pocket engagement
- compound 09: 0/100 orthosteric; 0/100 partial-pocket engagement
- compound 11: 10/100 orthosteric; qualifying pattern in 5/5 seeds

Selected compound 11 pose:
- run01, mode06
- Vinardo score: -3.5860 kcal/mol
- centroid distance from OHN: 4.7025 A

LasR PLIP:
- PLIP run records detected: 1
- Ala127 hydrophobic contact: 3.68 A
- Asp65 hydrogen-bond assignment: H-A 3.12 A; D-A 4.02 A; angle 148.51 degrees
- Asp65 contact retained as marginal/long, not a strong anchor

## PqsR exploratory docking

- compounds: 08, 09, 11
- independent runs: 5 per compound
- run-record files detected: 15
- all results explicitly exploratory because redocking failed

Selected representatives:
- compound 08: Vinardo -2.2480 kcal/mol; QZN centroid distance 9.4680 A
- compound 09: Vinardo -4.1040 kcal/mol; QZN centroid distance 6.8773 A
- compound 11: Vinardo -4.8680 kcal/mol; QZN centroid distance 4.3958 A

PqsR PLIP:
- PLIP run records detected: 3
- compound 11: Leu197/Leu207 profile; no Tyr258-centered interaction
- compound 09: strongest credible contact to Glu151:B
- compound 08: broad interfacial contact network
- chemically invalid Lys266 'carboxylate salt bridges' for neutral compounds 8/9 excluded

## Visual QC

- PqsR compounds 08, 09, and 11: supervised ChimeraX QC completed and closed
- LasR compound 11: supervised ChimeraX QC completed and closed
- no positive LasR panels generated for compounds 08 or 09
- receptor-specific visual-QC closure manifests verified before global closure

## Dedicated manifests created at global closure

- `99_logs_versions_checksums/checksums/sha256_step9_PqsR_exploratory_PLIP_2026-07-19.txt`
- `99_logs_versions_checksums/checksums/sha256_step9_redocking_complete_2026-07-19.txt`
- `99_logs_versions_checksums/checksums/sha256_step9_GLOBAL_CLOSURE_2026-07-19.txt`

## Interpretation boundary

Docking and PLIP provide structural hypotheses. They do not establish experimental
affinity, receptor occupancy, intracellular exposure, or causal inhibition of quorum sensing.
PqsR results remain exploratory.
