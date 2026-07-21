#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

ROOT=/data/BIOFILM_INSILICO_SUPPORT
CFGDIR="$ROOT/03_docking/vina_configs/production_LasR_2026-07-16"
OUTROOT="$ROOT/03_docking/vina_outputs/production/LasR_2UV0_2026-07-16"
MASTER="$OUTROOT/LasR_2UV0_production_master_summary.tsv"

printf 'configuration\tligand\tseed\tstart\tend\texit_status\tmode1_score_kcal_mol\telapsed\tmax_rss_kb\tpose_sha256\tlog_sha256\ttime_sha256\trecord_sha256\n' > "$MASTER"

RUNCOUNT=0

while IFS= read -r CONF; do
    NAME=$(basename "$CONF" .conf)
    LIGAND=$(awk -F' = ' '$1=="ligand"{print $2}' "$CONF")
    SEED=$(awk -F' = ' '$1=="seed"{print $2}' "$CONF")
    OUT=$(awk -F' = ' '$1=="out"{print $2}' "$CONF")
    DIR=$(dirname "$OUT")
    STEM="${OUT%_poses.pdbqt}"
    LOG="${STEM}_run.log"
    TIMELOG="${STEM}_time.log"
    RECORD="${STEM}_run_record.txt"

    test -s "$CONF"
    test -s "$LIGAND"
    test -d "$DIR"

    for F in "$OUT" "$LOG" "$TIMELOG" "$RECORD"; do
        test ! -e "$F"
    done

    START=$(date --iso-8601=seconds)

    {
        printf 'STEP 9 — LasR production docking\n'
        printf 'Configuration: %s\n' "$CONF"
        printf 'Ligand: %s\n' "$LIGAND"
        printf 'Seed: %s\n' "$SEED"
        printf 'Start: %s\n' "$START"
        printf 'Host: %s\n' "$(hostname)"
        printf 'Command: /usr/bin/vina --config %s\n' "$CONF"
        printf 'Vina version: '
        /usr/bin/vina --version 2>&1 | head -n 1
        printf '\nConfiguration SHA-256:\n'
        sha256sum "$CONF"
        printf '\nInput SHA-256:\n'
        sha256sum \
          "$ROOT/03_docking/receptors_pdbqt/LasR_2UV0.pdbqt" \
          "$LIGAND"
    } > "$RECORD"

    printf '\n===== STARTING %s =====\n' "$NAME"

    if /usr/bin/time -v -o "$TIMELOG" \
        /usr/bin/vina --config "$CONF" > "$LOG" 2>&1; then
        STATUS=0
    else
        STATUS=$?
    fi

    END=$(date --iso-8601=seconds)

    {
        printf '\nEnd: %s\n' "$END"
        printf 'Exit status: %s\n' "$STATUS"
        if [ -e "$OUT" ]; then
            printf '\nOutput SHA-256:\n'
            sha256sum "$OUT" "$LOG" "$TIMELOG"
        fi
    } >> "$RECORD"

    if [ "$STATUS" -ne 0 ]; then
        printf '\nERROR IN %s — LAST 50 LOG LINES:\n' "$NAME"
        tail -n 50 "$LOG"
        exit "$STATUS"
    fi

    MODE1=$(awk '/^   1[[:space:]]/{print $2;exit}' "$LOG")
    ELAPSED=$(awk -F': ' '/Elapsed \(wall clock\)/{print $2}' "$TIMELOG")
    MAXRSS=$(awk -F': ' '/Maximum resident set size/{print $2}' "$TIMELOG")

    POSE_SHA=$(sha256sum "$OUT" | awk '{print $1}')
    LOG_SHA=$(sha256sum "$LOG" | awk '{print $1}')
    TIME_SHA=$(sha256sum "$TIMELOG" | awk '{print $1}')
    RECORD_SHA=$(sha256sum "$RECORD" | awk '{print $1}')

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$NAME" "$(basename "$LIGAND")" "$SEED" "$START" "$END" "$STATUS" \
      "$MODE1" "$ELAPSED" "$MAXRSS" "$POSE_SHA" "$LOG_SHA" "$TIME_SHA" \
      "$RECORD_SHA" >> "$MASTER"

    RUNCOUNT=$((RUNCOUNT + 1))
    printf 'COMPLETED %d/15  %s  mode1=%s kcal/mol  elapsed=%s\n' \
      "$RUNCOUNT" "$NAME" "$MODE1" "$ELAPSED"

done < <(find "$CFGDIR" -maxdepth 1 -type f -name '*.conf' | sort)

test "$RUNCOUNT" -eq 15

printf '\nALL LASR PRODUCTION RUNS COMPLETED: %d/15\n' "$RUNCOUNT"
printf '\nMASTER SUMMARY SHA-256:\n'
sha256sum "$MASTER"
