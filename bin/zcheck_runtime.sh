#!/bin/bash

# zcheck_runtime.sh: Estimate total runtime of AMS MD simulations
# Usage: zcheck_runtime.sh [--table FILE] [path]
#   path: optional directory or log file. If a directory is given, scans recursively for
#         ams.results/ams.log files and reports each simulation plus a summary table.
#         If omitted, defaults to current directory.

set -euo pipefail

table_output=""

while [[ ${1:-} == --table ]]; do
  case "$1" in
    --table)
      shift
      table_output="${1-}"
      if [[ -z "$table_output" ]]; then
        echo "Error: --table requires a file path" >&2
        exit 2
      fi
      shift
      ;;
  esac
done

format_duration() {
  local total=${1:-0}
  local days=$(( total/86400 ))
  local hours=$(( (total%86400)/3600 ))
  local minutes=$(( (total%3600)/60 ))
  local seconds=$(( total%60 ))
  if (( days > 0 )); then
    printf "%dd %02d:%02d:%02d" "$days" "$hours" "$minutes" "$seconds"
  else
    printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds"
  fi
}

arg="${1-}"
declare -a logfiles=()

if [[ -z "$arg" ]]; then
  search_root="."
elif [[ -f "$arg" ]]; then
  case "$arg" in
    */ams.results/ams.log) logfiles=("$arg") ;;
    *) echo "Error: file must be an ams.results/ams.log" >&2; exit 1 ;;
  esac
elif [[ -d "$arg" ]]; then
  search_root="$arg"
else
  echo "Error: path not found: $arg" >&2
  exit 1
fi

if [[ ${#logfiles[@]} -eq 0 ]]; then
  mapfile -t logfiles < <(find "$search_root" -type f -path "*/ams.results/ams.log" | sort)
fi

if [[ ${#logfiles[@]} -eq 0 ]]; then
  echo "No ams.results/ams.log files found under ${search_root:-$arg}" >&2
  exit 1
fi

declare -A percent_table
declare -A elapsed_table
declare -A row_seen
declare -A col_seen
rows=()
cols=()

process_log() {
  local logfile="$1"

  echo "Using log file: $logfile"

  if [[ ! -f "$logfile" ]]; then
    echo "File not found: $logfile" >&2
    return 1
  fi

  local start_line end_line
  start_line=$(grep -E ">\s+0\s+0\.00\s+" "$logfile" | head -n1 || true)
  end_line=$(grep -E ">\s+[0-9]+\s+[0-9]+\.[0-9]+\s+" "$logfile" | tail -n1 || true)

  if [[ -z "$start_line" ]]; then
    echo "Error: Could not find simulation start line in log file" >&2
    return 1
  fi

  if [[ -z "$end_line" ]]; then
    echo "Error: Could not find any simulation step lines in log file" >&2
    return 1
  fi

  local raw_start_date start_time raw_end_date end_time
  raw_start_date=$(awk '{print $1}' <<<"$start_line" | tr -d '<>')
  start_time=$(awk '{print $2}' <<<"$start_line" | tr -d '<>')
  raw_end_date=$(awk '{print $1}' <<<"$end_line" | tr -d '<>')
  end_time=$(awk '{print $2}' <<<"$end_line" | tr -d '<>')

  local start_date end_date
  start_date=$(sed -E 's/^([A-Za-z]+)([0-9]+)-([0-9]+)/\1 \2 \3/' <<<"$raw_start_date")
  end_date=$(sed -E 's/^([A-Za-z]+)([0-9]+)-([0-9]+)/\1 \2 \3/' <<<"$raw_end_date")

  local start_epoch end_epoch elapsed
  start_epoch=$(date -d "${start_date} ${start_time}" +%s)
  end_epoch=$(date -d "${end_date} ${end_time}" +%s)
  elapsed=$(( end_epoch - start_epoch ))

  if (( elapsed < 0 )); then
    echo "Error: negative elapsed time computed for $logfile" >&2
    return 1
  fi

  local days hours minutes seconds
  days=$(( elapsed/86400 ))
  hours=$(( (elapsed%86400)/3600 ))
  minutes=$(( (elapsed%3600)/60 ))
  seconds=$(( elapsed%60 ))

  local sample0_time sample1_line sample1_time samplesize
  sample0_time=$(awk '{print $4}' <<<"$start_line")
  sample1_line=$(grep -E ">\s+[0-9]+\s+[0-9]+\.[0-9]+\s+" "$logfile" | sed -n '2p' || true)
  if [[ -z "$sample1_line" ]]; then
    echo "Error: Not enough MD steps to determine sampling frequency" >&2
    return 1
  fi

  sample1_time=$(awk '{print $4}' <<<"$sample1_line")
  samplesize=$(awk -v a="$sample1_time" -v b="$sample0_time" 'BEGIN{printf "%.2f", a-b}')

  if [[ $(echo "$samplesize <= 0" | bc -l) -eq 1 ]]; then
    echo "Error: Invalid sample size ($samplesize fs). Cannot determine time step." >&2
    return 1
  fi

  local last_time nsamples
  last_time=$(awk '{print $4}' <<<"$end_line")
  nsamples=$(awk -v lt="$last_time" -v ss="$samplesize" 'BEGIN{if (ss==0) {print 0} else {printf "%.2f", lt/ss}}')

  if [[ $(echo "$nsamples <= 0" | bc -l) -eq 1 ]]; then
    echo "Error: Invalid number of samples ($nsamples). Cannot calculate runtime." >&2
    return 1
  fi

  local time_per_sample
  time_per_sample=$(awk -v e="$elapsed" -v n="$nsamples" 'BEGIN{if (n==0){print 0} else {printf "%.2f", e/n}}')

  local target=1000000 needed_samples predicted_time total_seconds
  needed_samples=$(awk -v t="$target" -v ss="$samplesize" 'BEGIN{if (ss==0){print 0} else {printf "%.2f", t/ss}}')
  predicted_time=$(awk -v ns="$needed_samples" -v tps="$time_per_sample" 'BEGIN{printf "%.2f", ns*tps}')
  total_seconds=$(printf "%.0f" "$predicted_time")

  local pred_days pred_hours pred_minutes pred_seconds
  pred_days=$(( total_seconds/86400 ))
  pred_hours=$(( (total_seconds%86400)/3600 ))
  pred_minutes=$(( (total_seconds%3600)/60 ))
  pred_seconds=$(( total_seconds%60 ))

  local percent_done
  percent_done=$(awk -v ns="$nsamples" -v need="$needed_samples" 'BEGIN{if (need==0){print "0.00"} else {printf "%.2f", (ns/need)*100}}')

  echo "----------------------------------------"
  echo "Start timestamp: ${start_date} ${start_time}"
  echo "End timestamp:   ${end_date} ${end_time}"
  echo "Elapsed real time: $elapsed seconds"
  echo "Sample size: ${samplesize} fs"
  echo "Samples performed: $nsamples"
  echo "Avg. time per sample: $time_per_sample seconds"
  echo "Needed samples: $needed_samples samples"
  echo "----------------------------------------"
  printf "Current time taken for %s %% of the simulation: %d seconds (~%02d-%02d:%02d:%02d)\n" \
    "$percent_done" $elapsed $days $hours $minutes $seconds
  printf "Predicted runtime for %d fs: %d seconds (~%02d-%02d:%02d:%02d)\n" \
    $target $total_seconds $pred_days $pred_hours $pred_minutes $pred_seconds
  echo "----------------------------------------"

  local run_dir sim_dir sim_row run_name temp freq column
  run_dir=$(dirname "$logfile")          # .../ams.results
  run_dir=$(dirname "$run_dir")          # .../MD_xxx_300K_100fs
  sim_dir=$(dirname "$run_dir")          # .../MD_xxx
  sim_row=$(basename "$sim_dir")
  run_name=$(basename "$run_dir")
  temp=$(grep -oE '[0-9]+K(_[0-9]+)?' <<<"$run_name" | head -n1)
  freq=$(grep -oE '[0-9]+fs' <<<"$run_name" | tail -n1)
  column="${temp}_${freq:-unknown}"

  if [[ -n "$sim_row" && -z "${row_seen[$sim_row]:-}" ]]; then
    rows+=("$sim_row")
    row_seen["$sim_row"]=1
  fi

  if [[ -n "$column" && -z "${col_seen[$column]:-}" ]]; then
    cols+=("$column")
    col_seen["$column"]=1
  fi

  local key="${sim_row}|${column}"
  percent_table["$key"]="$percent_done"
  elapsed_table["$key"]="$elapsed"
}

for logfile in "${logfiles[@]}"; do
  process_log "$logfile" || true
done

if [[ ${#percent_table[@]} -eq 0 ]]; then
  exit 0
fi

readarray -t sorted_rows < <(printf '%s\n' "${rows[@]}" | sort)
readarray -t sorted_cols < <(printf '%s\n' "${cols[@]}" | sort -t_ -k1,1V -k2,2V)

table_buffer=$'\nSummary Table (Percent Complete / Elapsed Real Time)\n'

table_buffer+=$(printf '%-25s' "Simulation")
for col in "${sorted_cols[@]}"; do
  table_buffer+=$(printf ' | %-28s' "$col")
done
table_buffer+=$'\n'
table_buffer+=$(printf '%0.s-' {1..25})
for _ in "${sorted_cols[@]}"; do
  table_buffer+="+"
  table_buffer+=$(printf '%0.s-' {1..29})
done
table_buffer+=$'\n'

for row in "${sorted_rows[@]}"; do
  table_buffer+=$(printf '%-25s' "$row")
  for col in "${sorted_cols[@]}"; do
    key="${row}|${col}"
    if [[ -n "${percent_table[$key]:-}" ]]; then
      percent_val=${percent_table[$key]}
      elapsed_val=${elapsed_table[$key]}
      elapsed_fmt=$(format_duration "${elapsed_val}")
  cell=$(printf '%6.2f%% / %s' "$percent_val" "$elapsed_fmt")
      table_buffer+=$(printf ' | %-28s' "$cell")
    else
      table_buffer+=$(printf ' | %-28s' "--")
    fi
  done
  table_buffer+=$'\n'
done

printf '%s' "$table_buffer"

if [[ -n "$table_output" ]]; then
  printf "%s" "$table_buffer" >"$table_output"
  echo "Summary table written to $table_output"
fi
