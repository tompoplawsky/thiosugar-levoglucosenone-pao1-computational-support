# Documentation correction: QM geometry and CHELPG provenance

- Date: 2026-07-16
- Applied at: 2026-07-16T10:57:36+02:00
- Project root: `/data/BIOFILM_INSILICO_SUPPORT`

## Scope

Updated three living project documents after completion and QC of STEPS 3–4.
No raw ORCA outputs, geometries, CHELPG outputs, ADMET data, receptor files, or docking files were modified.

## Scientific corrections

- `compound_08`: GFN2-xTB was a pre-optimization step only; the final geometry is a converged omegaB97X-D3/def2-SVP/CPCM(water) minimum after a Cartesian-Hessian restart, with 0 imaginary frequencies.
- All three current frontier-orbital datasets are at the same finalized DFT level and are directly comparable within this dataset.
- `compound_11`: the historical CHELPG N-4 value of approximately -1.087 e was not reproduced and is not a current result or QC target.
- Current validated `compound_11` mapping: N-4 +0.31321861 e; N-3 -0.45537399 e; exocyclic amino N -0.95315983 e.

## Versioning and manifests

- Exact pre-update copies were archived under `99_logs_versions_checksums/documentation_snapshots/2026-07-16_pre_step8_qm_chelpg_update`.
- The historical STEP 4 checksum manifest was intentionally not modified; it remains a point-in-time record.
- Current corrected documents and archived snapshots are covered by `99_logs_versions_checksums/checksums/sha256_documentation_correction_2026-07-16.txt`.

## SHA-256 changes

| File | Pre-update SHA-256 | Corrected SHA-256 |
|---|---|---|
| `Ai_Instruction_PROJECT_UPDATED.md` | `f35b373a37f3550e085e2a650d12f7e0e911b9f9a08ca9d2b6c12e204336b7f4` | `4a7111102baafa5e781cf9ea02b55ee639e48bd584db912ae2adf5a2cf6c0752` |
| `Srodowisko_obliczeniowe_PUBLIC_REPO.md` | `5483ddffd43a1c215feea0b65afc0bf7b4ca5695195fca9a7a7e6d7cc4fe53c3` | `33d068c4ba12f6d49a5eed911d4defa36722dfbb3ffa10d69b306a5ddaf0eaed` |
| `PIPELINE.md` | `65cb93f38e138a18c96536f7a4ab6da13c2db2694424fabf35396c5fc2e26d67` | `e9ad38389002d58069325335218eee72502cfbb1228f21df54bc3ef563254849` |
