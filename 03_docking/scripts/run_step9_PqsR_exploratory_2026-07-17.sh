#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

ROOT="/data/BIOFILM_INSILICO_SUPPORT"
CFGDIR="$ROOT/03_docking/vina_configs/exploratory_PqsR_2026-07-17"
OUTROOT="$ROOT/03_docking/vina_outputs/exploratory/PqsR_4JVI_2026-07-17"
MASTER_LOG="$OUTROOT/PqsR_4JVI_exploratory_master_runner.log"
SUMMARY="$OUTROOT/PqsR_4JVI_exploratory_master_summary.tsv"
PROTOCOL="$ROOT/03_docking/vina_configs/step9_production_and_exploratory_docking_protocol_preregistered_2026-07-16.txt"
ADDENDUM="$ROOT/03_docking/vina_configs/step9_PqsR_exploratory_analysis_addendum_preregistered_2026-07-17.txt"
VINA="/usr/bin/vina"

test -x "$VINA"
test -d "$CFGDIR"
test -d "$OUTROOT"
test -s "$PROTOCOL"
test -s "$ADDENDUM"
test ! -e "$MASTER_LOG"
test ! -e "$SUMMARY"

exec > >(tee -a "$MASTER_LOG") 2>&1

printf 'configuration\tligand\tseed\tstart\tend\texit_status\tmode1_score_kcal_mol\telapsed\tmax_rss_kb\tpose_sha256\tlog_sha256\ttime_sha256\trecord_sha256\n' > "$SUMMARY"

echo "STEP 9 — PqsR 4JVI EXPLORATORY DOCKING"
echo "Started: $(date --iso-8601=seconds)"
echo "Interpretation status: exploratory; formal receptor-specific validation not achieved"
echo "Primary protocol SHA-256: $(sha256sum "$PROTOCOL" | awk '{print $1}')"
echo "Exploratory addendum SHA-256: $(sha256sum "$ADDENDUM" | awk '{print $1}')"
echo "Vina: $("$VINA" --version 2>&1 | tr '\n' ' ')"
echo

entries=(
  "compound_11 01 2026071691"
  "compound_11 02 2026071692"
  "compound_11 03 2026071693"
  "compound_11 04 2026071694"
  "compound_11 05 2026071695"
  "compound_09 01 2026071681"
  "compound_09 02 2026071682"
  "compound_09 03 2026071683"
  "compound_09 04 2026071684"
  "compound_09 05 2026071685"
  "compound_08 01 2026071671"
  "compound_08 02 2026071672"
  "compound_08 03 2026071673"
  "compound_08 04 2026071674"
  "compound_08 05 2026071675"
)

completed=0
total="${#entries[@]}"

for entry in "${entries[@]}"; do
    read -r ligand run seed <<< "$entry"

    name="${ligand}_PqsR_4JVI_run${run}_seed${seed}"
    config="$CFGDIR/${name}.conf"
    rundir="$OUTROOT/$ligand/run${run}_seed${seed}"
    pose="$rundir/${name}_poses.pdbqt"
    runlog="$rundir/${name}_run.log"
    timelog="$rundir/${name}_time.log"
    record="$rundir/${name}_run_record.txt"

    test -s "$config"
    test ! -e "$rundir"
    mkdir -p "$rundir"

    start="$(date --iso-8601=seconds)"

    echo "STARTING $((completed + 1))/$total: $name"
    echo "Start: $start"

    set +e
    /usr/bin/time -v -o "$timelog" \
        "$VINA" --config "$config" --out "$pose" \
        > "$runlog" 2>&1
    status=$?
    set -e

    end="$(date --iso-8601=seconds)"

    if [[ "$status" -ne 0 ]]; then
        {
            echo "PqsR exploratory docking run record"
            echo "Status: FAILED"
            echo "Configuration: $name"
            echo "Start: $start"
            echo "End: $end"
            echo "Exit status: $status"
            echo "Config: $config"
            echo "Run log: $runlog"
            echo "Time log: $timelog"
        } > "$record"

        echo "FAILED: $name, exit status $status"
        exit "$status"
    fi

    test -s "$pose"
    test -s "$runlog"
    test -s "$timelog"

    mode1="$(awk '/^REMARK VINA RESULT:/{print $4; exit}' "$pose")"
    elapsed="$(sed -n 's/^[[:space:]]*Elapsed (wall clock) time (h:mm:ss or m:ss):[[:space:]]*//p' "$timelog")"
    maxrss="$(sed -n 's/^[[:space:]]*Maximum resident set size (kbytes):[[:space:]]*//p' "$timelog")"

    test -n "$mode1"
    test -n "$elapsed"
    test -n "$maxrss"

    pose_sha="$(sha256sum "$pose" | awk '{print $1}')"
    log_sha="$(sha256sum "$runlog" | awk '{print $1}')"
    time_sha="$(sha256sum "$timelog" | awk '{print $1}')"

    {
        echo "PqsR exploratory docking run record"
        echo
        echo "Interpretation status: exploratory only"
        echo "Formal PqsR redocking validation: not achieved"
        echo "Configuration: $name"
        echo "Ligand: $ligand"
        echo "Receptor: PqsR_4JVI"
        echo "Seed: $seed"
        echo "Start: $start"
        echo "End: $end"
        echo "Exit status: $status"
        echo "Mode-1 Vinardo score: $mode1 kcal/mol"
        echo "Elapsed: $elapsed"
        echo "Maximum RSS: $maxrss kB"
        echo
        echo "Config: $config"
        echo "Pose output: $pose"
        echo "Run log: $runlog"
        echo "Time log: $timelog"
        echo
        echo "Protocol SHA-256: $(sha256sum "$PROTOCOL" | awk '{print $1}')"
        echo "Addendum SHA-256: $(sha256sum "$ADDENDUM" | awk '{print $1}')"
        echo "Config SHA-256: $(sha256sum "$config" | awk '{print $1}')"
        echo "Pose SHA-256: $pose_sha"
        echo "Run-log SHA-256: $log_sha"
        echo "Time-log SHA-256: $time_sha"
        echo
        echo "Scientific restriction:"
        echo "This run is hypothesis-generating and must not be described as"
        echo "validated PqsR binding, binding affinity, or confirmed target engagement."
    } > "$record"

    record_sha="$(sha256sum "$record" | awk '{print $1}')"

    printf '%s\t%s.pdbqt\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$name" "$ligand" "$seed" "$start" "$end" "$status" \
        "$mode1" "$elapsed" "$maxrss" \
        "$pose_sha" "$log_sha" "$time_sha" "$record_sha" \
        >> "$SUMMARY"

    completed=$((completed + 1))
    echo "COMPLETED $completed/$total: $name"
    echo
done

echo "ALL PQSR EXPLORATORY RUNS COMPLETED: $completed/$total"
echo
echo "MASTER SUMMARY SHA-256:"
sha256sum "$SUMMARY"
