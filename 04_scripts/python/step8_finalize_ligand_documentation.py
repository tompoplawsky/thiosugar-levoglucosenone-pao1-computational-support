#!/usr/bin/env python3
"""
Finalize STEP 8 reproducibility documentation.

Default: dry-run only.
Apply: copy the STEP 8 preparation/audit scripts into 04_scripts/python,
record their paths and hashes in the protocol and run log, and regenerate
the STEP 8 checksum manifest.

The script refuses ambiguous or previously finalized states.
It does not update backups.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path("/data/BIOFILM_INSILICO_SUPPORT")
DESKTOP = Path.home() / "Desktop"

PROTOCOL = ROOT / "03_docking/ligand_preparation/step8_ligand_preparation_protocol.json"
RUN_LOG = ROOT / "99_logs_versions_checksums/run_logs/step8_ligand_pdbqt_preparation_2026-07-16.md"
VERSION_FILE = ROOT / "99_logs_versions_checksums/software_versions/step8_software_versions_2026-07-16.txt"
MANIFEST = ROOT / "99_logs_versions_checksums/checksums/sha256_step8_ligands_2026-07-16.txt"

EXPECTED_PRE_HASHES = {
    PROTOCOL: "bd22473ae6178afe70f2e3c04b544a172a8200b443a9aa933615d6cd35cff13d",
    RUN_LOG: "eda06ff596216ddc4bc3417d3b8f331eb6ec90f8186462b6b23b373e8d9f0f01",
    VERSION_FILE: "b750b9c7caa9089c424673d4e58e59c5095d4d89071d7c2d87276da598bb86f0",
    MANIFEST: "e6f38ad471dcd97ebd09ea056d1f3743f318330c9a2db161855e6850c46237cd",
}

SCRIPT_NAMES = [
    "step8b_audit_ligand_topology_mapping.py",
    "step8b_map_meeko_rotors.py",
    "step8b_validate_proposed_torsion_rules.py",
    "step8_prepare_final_qm_ligands.py",
    "step8_finalize_ligand_documentation.py",
]

DEST_SCRIPT_DIR = ROOT / "04_scripts/python"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_manifest_entries(path: Path) -> list[Path]:
    entries: list[Path] = []
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        try:
            digest, filename = line.split("  ", 1)
        except ValueError as exc:
            raise RuntimeError(f"Malformed manifest line {lineno}: {line}") from exc
        file_path = Path(filename)
        if not file_path.is_absolute():
            raise RuntimeError(f"Manifest entry is not absolute: {file_path}")
        if not file_path.is_file():
            raise RuntimeError(f"Manifest entry missing: {file_path}")
        actual = sha256(file_path)
        if actual != digest:
            raise RuntimeError(
                f"Manifest validation failed for {file_path}: expected {digest}, got {actual}"
            )
        entries.append(file_path)
    if len(entries) != 12:
        raise RuntimeError(f"Expected 12 current manifest entries, found {len(entries)}")
    return entries


def audit_pre_state() -> tuple[list[Path], dict[str, dict[str, str]]]:
    for path, expected in EXPECTED_PRE_HASHES.items():
        if not path.is_file():
            raise RuntimeError(f"Missing required file: {path}")
        actual = sha256(path)
        if actual != expected:
            raise RuntimeError(
                f"Unexpected pre-finalization hash for {path}\n"
                f"Expected: {expected}\nActual:   {actual}"
            )

    manifest_entries = read_manifest_entries(MANIFEST)

    scripts: dict[str, dict[str, str]] = {}
    for name in SCRIPT_NAMES:
        source = DESKTOP / name
        destination = DEST_SCRIPT_DIR / name
        if not source.is_file():
            raise RuntimeError(f"Missing Desktop source script: {source}")
        if destination.exists():
            raise RuntimeError(f"Refusing to overwrite existing repository script: {destination}")
        scripts[name] = {
            "source": str(source),
            "destination": str(destination),
            "relative_destination": str(destination.relative_to(ROOT)),
            "sha256": sha256(source),
        }

    protocol_data = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if "reproducibility_scripts" in protocol_data:
        raise RuntimeError("Protocol already contains reproducibility_scripts")

    run_log_text = RUN_LOG.read_text(encoding="utf-8")
    if "## Reproducibility scripts" in run_log_text:
        raise RuntimeError("Run log already contains reproducibility scripts section")

    return manifest_entries, scripts


def updated_protocol(scripts: dict[str, dict[str, str]]) -> str:
    data = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    data["reproducibility_scripts"] = {
        name: {
            "path": item["relative_destination"],
            "sha256": item["sha256"],
            "role": {
                "step8b_audit_ligand_topology_mapping.py":
                    "Read-only audit of final XYZ atom order versus audited MOL/SDF topology.",
                "step8b_map_meeko_rotors.py":
                    "Read-only mapping of default Meeko branches to original ligand bonds.",
                "step8b_validate_proposed_torsion_rules.py":
                    "Read-only validation that chemistry-aware SMARTS rules remove only intended conjugated rotors.",
                "step8_prepare_final_qm_ligands.py":
                    "Canonical dry-run/apply preparation of final-QM SDF and ligand PDBQT files.",
                "step8_finalize_ligand_documentation.py":
                    "Finalization of script provenance and STEP 8 checksum manifest.",
            }[name],
        }
        for name, item in scripts.items()
    }
    return json.dumps(data, indent=2, ensure_ascii=False) + "\n"


def updated_run_log(scripts: dict[str, dict[str, str]]) -> str:
    text = RUN_LOG.read_text(encoding="utf-8").rstrip() + "\n\n"
    lines = [
        "## Reproducibility scripts",
        "",
        "| Script | Role | SHA-256 |",
        "|---|---|---|",
    ]
    roles = {
        "step8b_audit_ligand_topology_mapping.py":
            "Topology and atom-order audit",
        "step8b_map_meeko_rotors.py":
            "Default Meeko rotor mapping",
        "step8b_validate_proposed_torsion_rules.py":
            "Validation of chemistry-aware torsion rules",
        "step8_prepare_final_qm_ligands.py":
            "Canonical ligand SDF/PDBQT preparation",
        "step8_finalize_ligand_documentation.py":
            "Documentation and checksum finalization",
    }
    for name, item in scripts.items():
        lines.append(
            f"| `{item['relative_destination']}` | {roles[name]} | `{item['sha256']}` |"
        )
    lines += [
        "",
        "The preparation and audit scripts were archived in the repository before STEP 8 backup.",
        "The STEP 8 checksum manifest was regenerated after adding these scripts and updating this documentation.",
        "",
    ]
    return text + "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    manifest_entries, scripts = audit_pre_state()

    protocol_text = updated_protocol(scripts)
    run_log_text = updated_run_log(scripts)

    print("PRE-FINALIZATION AUDIT PASSED")
    print(f"Validated current manifest entries: {len(manifest_entries)}")
    print("Scripts to archive:")
    for name, item in scripts.items():
        print(f"  {item['sha256']}  {item['relative_destination']}")

    if not args.apply:
        print("\nDRY-RUN PASSED")
        print("Repository was not modified.")
        print("Planned documentation updates:")
        print(f"  - {PROTOCOL}")
        print(f"  - {RUN_LOG}")
        print(f"  - {MANIFEST}")
        return 0

    DEST_SCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    for name, item in scripts.items():
        source = Path(item["source"])
        destination = Path(item["destination"])
        shutil.copy2(source, destination)

    PROTOCOL.write_text(protocol_text, encoding="utf-8")
    RUN_LOG.write_text(run_log_text, encoding="utf-8")

    final_entries = list(manifest_entries)
    for name in SCRIPT_NAMES:
        final_entries.append(DEST_SCRIPT_DIR / name)

    # Ensure updated protocol/run log are represented only once.
    deduped: list[Path] = []
    seen: set[Path] = set()
    for path in final_entries:
        if path not in seen:
            deduped.append(path)
            seen.add(path)

    lines = [f"{sha256(path)}  {path}" for path in deduped]
    MANIFEST.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Final internal validation.
    validated = read_manifest_entries_flexible(MANIFEST, expected_count=len(deduped))

    print("\nAPPLY PASSED")
    print(f"Final manifest entries: {len(validated)}")
    for path in [*deduped, MANIFEST]:
        print(f"{sha256(path)}  {path}")
    print("Backups have not yet been updated.")
    return 0


def read_manifest_entries_flexible(path: Path, expected_count: int) -> list[Path]:
    entries: list[Path] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, filename = line.split("  ", 1)
        file_path = Path(filename)
        if sha256(file_path) != digest:
            raise RuntimeError(f"Final manifest verification failed: {file_path}")
        entries.append(file_path)
    if len(entries) != expected_count:
        raise RuntimeError(
            f"Final manifest count mismatch: expected {expected_count}, got {len(entries)}"
        )
    return entries


if __name__ == "__main__":
    raise SystemExit(main())
