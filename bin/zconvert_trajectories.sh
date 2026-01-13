#!/bin/bash
# =====================================================================
# zconvert_trajectories.sh
# ---------------------------------------------------------------------
# Can be used to run conversions locally or via SLURM
# Writes matti .xyz or .pdb files from AMS rkf trajectories using kf_converter.py
#
# Author: L. Kuijken
# Last updated: 2025-11-13
# =====================================================================


# =====================================================================
# INITIAL SETUP
# =====================================================================
PY="$HOME/master_thesis_project/Tools/kf_converter.py"
LOGDIR="$HOME/master_thesis_project/Trajectories/Convert_logs"
MASTERLOG="$LOGDIR/master_summary.log"
DATETIME=$(date '+%Y-%m-%d %H:%M:%S')
TIME=$(date +%H:%M:%S)

mkdir -p "$LOGDIR"

# =====================================================================
# SHOW HELP MESSAGE
# =====================================================================
show_help() {
    echo "Usage: $0 BASE_FOLDER [-f format] [-m mode]"
    echo
    echo "Convert AMS trajectories using either local or SLURM execution."
    echo
    echo "Arguments:"
    echo "  BASE_FOLDER          Base directory containing AMS simulation results"
    echo
    echo "Options:"
    echo "  -r, --region         Region to process (default: all regions)"
    echo "  -s, --sampling       Output sampling interval (e.g. 100fs)"
    echo "  -t, --max_time       Time range to process (default: all time)"
    echo "  -f, --format         Output format: 'matti' or 'pdb' (default: matti)"
    echo "  -o, --out            Output directory (default: ~/master_thesis_project/Trajectories/<format>)"
    echo "  -w, --overwrite      Overwrite mode: prompt, skip, or force (default: prompt)"
    echo "  -m, --mode           Execution mode: 'local' or 'slurm' (default: test)"
    echo "  -N, --molecules      Comma-separated list of water molecule counts (default: 33,65,98,130)"
    echo "  -H, --hydronium      Comma-separated list of hydronium counts (default: 0,1,2,3,4)"
    echo "  -T, --temperatures   Comma-separated list of temperatures (default: 300,350,400)"
    echo "  -I, --input_sampling Input sampling to filter (e.g., 2fs, 100fs, 0.1ps) - empty means all"
    echo "  -h, --help           Show this message and exit"
    exit 0
}

# =====================================================================
# CHECK FOR HELP FLAG AND BASE FOLDER BEFORE BACKGROUND EXECUTION
# =====================================================================
for arg in "$@"; do
    if [ "$arg" = "-h" ] || [ "$arg" = "--help" ]; then
        show_help
        exit 0
    fi
done

# Quick check if BASE_FOLDER is provided (must be first non-flag argument)
HAS_BASE=0
for arg in "$@"; do
    if [[ ! "$arg" =~ ^- ]]; then
        HAS_BASE=1
        break
    fi
done

if [ "$HAS_BASE" -eq 0 ]; then
    echo "Error: No BASE_FOLDER provided."
    echo ""
    show_help
    exit 1
fi

# =====================================================================
# BACKGROUND EXECUTION
# =====================================================================
echo "" >> "$MASTERLOG"
echo "[$DATETIME]: Running zconvert_trajectories.sh > $MASTERLOG..."
find "$LOGDIR" -type f ! -name '*master*.log' -delete

if [ -z "$CONVERT_BACKGROUND" ]; then
    export CONVERT_BACKGROUND=1
    exec "$0" "$@" >> "$MASTERLOG" 2>&1 &
    exit $?
fi

# Default parameters
RUN_MODE="test"
FORMAT="matti"
N_FILTER="33,65,98,130"
H_FILTER="0,1,2,3,4"
T_FILTER="300,350,400"
I_FILTER=""

# =====================================================================
# UTILITY FUNCTIONS
# =====================================================================
parse_time_to_fs() {
    local time_str="$1"
    local value="${time_str//[^0-9.]/}"  # Extract number
    local unit="${time_str//[0-9.]/}"    # Extract unit
    
    case "${unit,,}" in  # Convert to lowercase
        fs)
            echo "$value"
            ;;
        ps)
            echo "$(awk "BEGIN {print int($value * 1000)}")"
            ;;
        ns)
            echo "$(awk "BEGIN {print int($value * 1000000)}")"
            ;;
        *)
            echo "Error: Invalid time unit '$unit' in '$time_str'" >&2
            return 1
            ;;
    esac
}

# =====================================================================
# PARSE ARGUMENTS FROM COMMAND LINE
# =====================================================================
PYTHON_ARGS=""
while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help) show_help ;;
        -r|--region) PYTHON_ARGS="$PYTHON_ARGS --region $2"; shift 2 ;;
        -s|--sampling) PYTHON_ARGS="$PYTHON_ARGS --output_sampling $2"; shift 2 ;;
        -S|--max_time) PYTHON_ARGS="$PYTHON_ARGS --max_time $2"; shift 2 ;;
        -f|--format) FORMAT="$2"; PYTHON_ARGS="$PYTHON_ARGS --format $2"; shift 2 ;;
        -o|--out) PYTHON_ARGS="$PYTHON_ARGS --out $2"; shift 2 ;;
        -w|--overwrite) PYTHON_ARGS="$PYTHON_ARGS --overwrite $2"; shift 2 ;;
        -m|--mode) RUN_MODE="$2"; shift 2 ;;
        -N|--molecules) N_FILTER="$2"; shift 2 ;;
        -H|--hydronium) H_FILTER="$2"; shift 2 ;;
        -T|--temperatures) T_FILTER="$2"; shift 2 ;;
        -I|--input_sampling) I_FILTER="$2"; shift 2 ;;
        *)
            if [ -z "$BASE" ]; then
                BASE="$1"
            else
                echo "[$TIME] Unknown argument: $1"
                exit 1
            fi
            shift ;;
    esac
done

# Validation
if [ -z "$BASE" ]; then
    echo "[$TIME] Error: no base folder specified."
    exit 1
fi
if [ ! -d "$BASE" ]; then
    echo "[$TIME] Error: directory not found: $BASE"
    exit 1
fi

# Auto-adjust overwrite for SLURM
if [ "$RUN_MODE" = "slurm" ]; then
    if [[ "$PYTHON_ARGS" != *"--overwrite"* ]]; then
        PYTHON_ARGS="$PYTHON_ARGS --overwrite skip"
    fi
fi

# =====================================================================
# FOLDER LIST GENERATION
# =====================================================================
generate_filtered_list() {
    FOLDERLIST="$LOGDIR/folderlist.txt"

    echo "[$TIME] Scanning: $BASE"
    
    # Convert sampling filter if provided
    I_FS=""
    if [ -n "$I_FILTER" ]; then
        I_FS=$(parse_time_to_fs "$I_FILTER")
        EXIT_CODE=$?
        
        if [ $EXIT_CODE -ne 0 ] || [ -z "$I_FS" ]; then
            echo "[$TIME] Error: Invalid sampling value '$I_FILTER'"
            exit 1
        fi
        echo "[$TIME] Filters -> N=[$N_FILTER], H=[$H_FILTER], T=[$T_FILTER], I=[${I_FS}fs]"
    else
        echo "[$TIME] Filters -> N=[$N_FILTER], H=[$H_FILTER], T=[$T_FILTER], I=[prefer 100fs]"
    fi

    # Convert comma-separated filters to arrays
    IFS=',' read -r -a T_ARR <<< "$T_FILTER"
    IFS=',' read -r -a N_ARR <<< "$N_FILTER"
    IFS=',' read -r -a H_ARR <<< "$H_FILTER"

    # Build a combined pattern list
    TMPFILE=$(mktemp)
    for T in "${T_ARR[@]}"; do
        for N in "${N_ARR[@]}"; do
            for H in "${H_ARR[@]}"; do
                N_H2O=$((N - H))
                if [ -n "$I_FS" ]; then
                    # Include sampling in pattern - full exact match
                    echo "MD_${N_H2O}H2O_${H}H3O_${T}K_${I_FS}fs\$" >> "$TMPFILE"
                else
                    # Match any sampling - use looser pattern
                    echo "MD_${N_H2O}H2O_${H}H3O_${T}K_[0-9]+fs" >> "$TMPFILE"
                fi
            done
        done
    done

    # Convert to single regex joined by |
    PATTERN=$(paste -sd'|' "$TMPFILE")
    rm -f "$TMPFILE"

    # Search and filter
    find "$BASE" -type f -path "*/ams.results/ams.rkf" | sed 's|/ams.results/ams.rkf$||' |
    grep -E "$PATTERN" | sort -u > "$LOGDIR/tmp_allfolders.txt"

    # If no input filter specified, deduplicate: prefer 100fs if both 20fs and 100fs exist
    if [ -z "$I_FS" ]; then
        echo "[$TIME] Deduplicating: preferring 100fs over other sampling frequencies..."
        
        awk -F'/' '
        {
            folder=$0
            base=$NF
            key=base
            # Remove sampling suffix to create grouping key
            sub(/_[0-9]+fs$/, "", key)
            
            # Prefer 100fs, otherwise keep first occurrence
            if (!(key in best) || (folder ~ /_100fs/ && best[key] !~ /_100fs/)) {
                best[key]=folder
            }
        }
        END { for (k in best) print best[k] }' "$LOGDIR/tmp_allfolders.txt" | sort > "$FOLDERLIST"
        
        rm -f "$LOGDIR/tmp_allfolders.txt"
    else
        # Input filter specified, no deduplication needed
        mv "$LOGDIR/tmp_allfolders.txt" "$FOLDERLIST"
    fi

    NUMFOLDERS=$(wc -l < "$FOLDERLIST")
    if [ "$NUMFOLDERS" -eq 0 ]; then
        echo "No matching simulations found. Exiting." 
        exit 0
    fi

    echo "[$TIME] Found $NUMFOLDERS matching simulations." 
}

# =====================================================================
# LOCAL EXECUTION
# =====================================================================
run_local() {
    JOBLOG="$LOGDIR/convert_local_$(date +%Y%m%d_%H%M%S).log"

    while read -r FOLDER; do
        [ -z "$FOLDER" ] && continue
        BASENAME=$(basename "$FOLDER")
        LOGFILE="$LOGDIR/${BASENAME}_convert.log"
        echo "[$TIME] Converting: $FOLDER" | tee -a "$JOBLOG"
        amspython "$PY" "$FOLDER" $PYTHON_ARGS >> "$LOGFILE" 2>&1
        pids+=($!)
    done < "$FOLDERLIST"
}
# =====================================================================
# SLURM EXECUTION
# =====================================================================
run_slurm() {
    # Pre-create logs
    while read -r FOLDER; do
        BASENAME=$(basename "$FOLDER")
        echo "[$TIME] Scheduled: $FOLDER" > "$LOGDIR/${BASENAME}_convert.log"
    done < "$FOLDERLIST"

    # Submit SLURM array
    JOBID=$(sbatch --wait --job-name="convert" \
           --partition=mech.student.q \
           --nodes=1 --ntasks=1 --cpus-per-task=1 \
           --mem=4G --time=2:00:00 \
           --array="1-$NUMFOLDERS%8" \
           --output=/dev/null \
           --wrap="
                TASK_ID=\$SLURM_ARRAY_TASK_ID
                FOLDER=\$(sed -n \"\${TASK_ID}p\" $FOLDERLIST)
                BASENAME=\$(basename \"\$FOLDER\")
                LOGFILE=\"$LOGDIR/\${BASENAME}_convert.log\"

                echo \"[\$(date +%H:%M:%S)] Running KFConverter for \$FOLDER\" >> \"\$LOGFILE\"
                echo \"Python args: $PYTHON_ARGS\" >> \"\$LOGFILE\"
                amspython \"$PY\" \"\$FOLDER\" $PYTHON_ARGS >> \"\$LOGFILE\" 2>&1
                STATUS=\$?
                if [ \$STATUS -eq 0 ]; then
                    echo \"[\$(date +%H:%M:%S)] Completed successfully.\" >> \"\$LOGFILE\"
                else
                    echo \"[\$(date +%H:%M:%S)] FAILED (exit \$STATUS).\" >> \"\$LOGFILE\"
                fi
                echo \"[\$(date +%H:%M:%S)] Job finished.\" >> \"\$LOGFILE\"
           " | awk '{print $4}')
}

# =====================================================================
# MAIN EXECUTION
# =====================================================================
generate_filtered_list
if [ "$RUN_MODE" = "local" ]; then
    echo "[$TIME] Starting LOCAL conversion ($NUMFOLDERS sims)..."
    run_local
    echo "[$TIME] LOCAL conversion completed."
elif [ "$RUN_MODE" = "slurm" ]; then
    echo "[$TIME] Starting SLURM conversion ($NUMFOLDERS sims)..." 
    run_slurm
    echo "[$TIME] SLURM conversion completed." 
else
    echo "[$TIME] Test run completed." 
    echo "No conversion executed. Use -m local or -m slurm to run conversions." 
fi