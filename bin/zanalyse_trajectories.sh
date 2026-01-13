#!/bin/bash
# =====================================================================
# zanalyse_trajectories.sh
# ---------------------------------------------------------------------
# Parallel trajectory analysis launcher for LionAnalysis.
# Processes multiple trajectories with configurable parallelism.
#
# Usage:
#   zanalyse_trajectories.sh BASE_FOLDER [options]
#   nohup zanalyse_trajectories.sh BASE_FOLDER [options] &
#
# Author: L. Kuijken
# Last updated: 2025-11-13
# =====================================================================


# =====================================================================
# INITIAL SETUP
# =====================================================================
PY="$HOME/master_thesis_project/Tools/trajectory.py"
LOGDIR="$HOME/master_thesis_project/LionAnalysis/Analysis_logs"
MASTERLOG="$LOGDIR/master_summary.log"

# Function to get current timestamp
timestamp() {
    date '+%H:%M:%S'
}

mkdir -p "$LOGDIR"

# =====================================================================
# SHOW HELP MESSAGE
# =====================================================================
show_help() {
    echo "Usage: $0 BASE_FOLDER [options]"
    echo
    echo "Parallel trajectory analysis launcher for LionAnalysis."
    echo
    echo "Arguments:"
    echo "  BASE_FOLDER             Directory containing trajectories or single trajectory file"
    echo
    echo "Analysis Options:"
    echo "  --atom_xyz INDEX        Track specific atom (single index, e.g., '700' for O atom)"
    echo "  --rdf_water             Enable RDF analysis for water"
    echo "  --rdf_mof               Enable RDF analysis for MOF"
    echo "  --msd_water             Enable water mean square displacement"
    echo "  --msd_hydronium         Enable hydronium mean square displacement"
    echo "  --h_bonds_water         Enable hydrogen bond analysis for water"
    echo "  --h_bonds_mof           Enable hydrogen bond analysis for MOF"
    echo "  --ptfel_water           Enable PTFEL (proton transfer free energy landscape) for water"
    echo "  --ptfel_mof             Enable PTFEL analysis for MOF"
    echo "  --res_times             Enable residence time analysis"
    echo "  --conductivity          Enable proton conductivity calculation"
    echo "  --debug                 Enable debug diagnostics"
    echo "  --trajectory_snapshot   Enable trajectory snapshot plotting"
    echo
    echo "Execution Options:"
    echo "  -w, --overwrite MODE    Overwrite mode: prompt, skip, or force (default: skip)"
    echo "  -s, --stride N          Analysis stride (frames per step, default: 1)"
    echo "  -t, --threads N         Threads per job (default: 8)"
    echo "  -j, --jobs N            Number of parallel jobs (default: 40)"
    echo
    echo "Other Options:"
    echo "  -h, --help           Show this help message"
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
# Execute in background with all output to logs
echo "" >> "$MASTERLOG"
echo "[$(date '+%Y-%m-%d %H:%M:%S')]: Running zanalyse_trajectories.sh > $MASTERLOG..."
find "$LOGDIR" -type f ! -name '*master*.log' -delete

if [ -z "$ANALYSE_BACKGROUND" ]; then
    export ANALYSE_BACKGROUND=1
    exec "$0" "$@" >> "$MASTERLOG" 2>&1 &
    exit $?
fi

# Default parameters
MAX_JOBS=40
MAX_THREADS=160


# =====================================================================
# PARSE ARGUMENTS FROM COMMAND LINE
# =====================================================================
PYTHON_ARGS=""
BASE=""

while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help) show_help ;;
        -w|--overwrite) PYTHON_ARGS="$PYTHON_ARGS --overwrite $2"; shift 2 ;;
        -s|--stride) PYTHON_ARGS="$PYTHON_ARGS --stride $2"; shift 2 ;;
        -t|--threads) PYTHON_ARGS="$PYTHON_ARGS --n_threads $2"; shift 2 ;;
        -j|--jobs) MAX_JOBS="$2"; shift 2 ;;
        --atom_xyz) PYTHON_ARGS="$PYTHON_ARGS --atom_xyz $2"; shift 2 ;;
        --rdf_water|--rdf_mof|--msd_water|--msd_hydronium|--h_bonds_water|--h_bonds_mof|--ptfel_water|--ptfel_mof|--res_times|--conductivity|--debug|--trajectory_snapshot)
            PYTHON_ARGS="$PYTHON_ARGS $1"; shift ;;
        *)
            if [ -z "$BASE" ]; then
                BASE="$1"
            else
                echo "[$(timestamp)] Unknown argument: $1"
                exit 1
            fi
            shift ;;
    esac
done


# =====================================================================
# VALIDATION
# =====================================================================
if [ -z "$BASE" ]; then
    echo "[$(timestamp)] Error: No BASE_FOLDER provided."
    exit 1
fi

if [ ! -d "$BASE" ] && [ ! -f "$BASE" ]; then
    echo "[$(timestamp)] Error: Base folder or file not found: $BASE"
    exit 1
fi


# =====================================================================
# TRAJECTORY DISCOVERY
# =====================================================================
JOBLOG="$LOGDIR/analyse_$(date +%Y%m%d_%H%M%S).log"

if [ -f "$BASE" ]; then
    # Single trajectory file provided
    FILES=("$BASE")
else
    # Search directory for .xyz trajectories
    readarray -t FILES < <(find "$BASE" -type f -name "*_traj.xyz" 2>/dev/null | sort -u)
fi

TOTAL=${#FILES[@]}

if [ "$TOTAL" -eq 0 ]; then
    echo "[$(timestamp)] No .xyz trajectories found in $BASE" | tee -a "$JOBLOG"
    exit 0
fi


# =====================================================================
# PARALLEL EXECUTION
# =====================================================================
echo "[$(timestamp)] Found $TOTAL trajectories in $BASE" | tee -a "$JOBLOG"
echo "[$(timestamp)] Command: $0 $BASE $PYTHON_ARGS"
echo "[$(timestamp)] Analysis run started with $MAX_JOBS parallel jobs"

COUNT=0
FAILED=0

for FILE in "${FILES[@]}"; do
    COUNT=$((COUNT + 1))
    NAME=$(basename "$FILE")
    LOGFILE="$LOGDIR/${NAME%.xyz}_analyse.log"

    echo "[$(timestamp)] [$COUNT/$TOTAL] Launching: $NAME" | tee -a "$JOBLOG"

    # Launch analysis in background
    (
        amspython "$PY" "$FILE" $PYTHON_ARGS >> "$LOGFILE" 2>&1
        STATUS=$?
        if [ $STATUS -eq 0 ]; then
            echo "[$(timestamp)] [$COUNT/$TOTAL] Completed: $NAME" | tee -a "$JOBLOG"
        else
            echo "[$(timestamp)] [$COUNT/$TOTAL] FAILED: $NAME (exit $STATUS)" | tee -a "$JOBLOG"
            echo "$NAME" >> "$LOGDIR/failed_analyses.txt"
        fi
    ) &

    # Job limit control
    while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do
        sleep 2
    done
done

# Wait for all background jobs to complete
wait


# =====================================================================
# SUMMARY
# =====================================================================
if [ -f "$LOGDIR/failed_analyses.txt" ]; then
    FAILED=$(wc -l < "$LOGDIR/failed_analyses.txt")
fi

SUCCESSFUL=$((TOTAL - FAILED))

echo "===================================================================="
echo "[$(timestamp)] Analysis run completed"
echo "[$(timestamp)] Total trajectories: $TOTAL"
echo "[$(timestamp)] Successful: $SUCCESSFUL"
echo "[$(timestamp)] Failed: $FAILED"
echo "===================================================================="

if [ "$FAILED" -gt 0 ]; then
    echo "[$(timestamp)] Failed analyses listed in: $LOGDIR/failed_analyses.txt"
fi