#!/bin/bash
# =====================================================================
# zrun_simulations.sh
# ---------------------------------------------------------------------
# Create and submit AMS MD simulations locally or via SLURM.
# Loops through parameter combinations to create individual simulation
# folders and optionally submit them with thread limit control.
#
# Usage:
#   zrun_simulations.sh [options]
#
# Author: L. Kuijken
# Last updated: 2025-11-13
# =====================================================================


# =====================================================================
# INITIAL SETUP
# =====================================================================
PROJECT_DIR="$HOME/master_thesis_project"
PY="$PROJECT_DIR/Tools/simulation.py"
LOGDIR="$PROJECT_DIR/AMS_simulations/Simulation_logs"
MASTERLOG="$LOGDIR/master_summary.log"

# Function to get current timestamp
timestamp() {
    date '+%H:%M:%S'
}

mkdir -p "$LOGDIR"
echo "" >> "$MASTERLOG"
echo "[$DATETIME]: Running zrun_simulations.sh > $MASTERLOG..."
find "$LOGDIR" -type f ! -name '*master*.log' -delete

# =====================================================================
# SHOW HELP MESSAGE
# =====================================================================
show_help() {
    echo "Usage: $0 -N NUM_MOLS [options]"
    echo
    echo "Create and optionally submit AMS MD simulations."
    echo
    echo "Required Arguments:"
    echo "  -N, --num_mols N     Total number of molecules (H2O + H3O)"
    echo
    echo "Options:"
    echo "  -H, --hydronium      Comma-separated H3O counts (default: 0,1,2,3,4)"
    echo "  -T, --temperatures   Comma-separated temperatures in K (default: 300,350,400)"
    echo "  -S, --n_steps N      Number of MD steps (default: 10000)"
    echo "  -f, --sampling N     Sampling frequency (default: 100)"
    echo "  -j, --threads N      Threads per simulation (default: 32)"
    echo "  -J, --max_threads N  Max total threads for MOPS/HPC (default: 128)"
    echo "  -m, --machine TYPE   Machine type: MOPS, HPC, or SNELLIUS (default: MOPS)"
    echo "  --START              Actually submit simulations (default: test mode)"
    echo "  --base_dir PATH      Base directory for simulations (default: $PROJECT_DIR/AMS_simulations)"
    echo "  --mol_path PATH      Path to molecule files (default: $PROJECT_DIR/AMS_simulations/molecules)"
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
if [ -z "$SIM_BACKGROUND" ]; then
    export SIM_BACKGROUND=1
    exec "$0" "$@" >> "$MASTERLOG" 2>&1 &
    exit $?
fi
 
# Default parameters
RUN_MODE="test"
MACHINE="MOPS"
N_STEPS=10000
SAMPLING_FREQ=100
SIM_THREADS=32
MAX_THREADS=128
T_FILTER="300,350,400"
H_FILTER="0,1,2,3,4"


# =====================================================================
# SUBMISSION FUNCTIONS
# =====================================================================
submit_mops() {
    local sim_folder="$1"
    local sim_name="$2"
    
    echo "[$(timestamp)] Submitting (local): $sim_name"
    cd "$sim_folder"
    export OMP_NUM_THREADS=$SIM_THREADS
    nohup ./run.run > out 2>&1 &
    echo $! > jobid_$!.jid
    pwd >> jobid_$!.jid
    date >> jobid_$!.jid
    echo "[$(timestamp)] $sim_folder $! [go]" >> "$HOME/log_adf"
    cd - > /dev/null
}

submit_hpc() {
    local sim_folder="$1"
    local sim_name="$2"
    
    echo "[$(timestamp)] Submitting (SLURM/HPC): $sim_name"
    sbatch \
        --job-name="$sim_name" \
        --partition=thin \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=$SIM_THREADS \
        --mem-per-cpu=2G \
        --time=5-00:00:00 \
        --output="$sim_folder/%x_%j.out" \
        --error="$sim_folder/%x_%j.err" \
        --export=ALL,SIM_FOLDER="$sim_folder",MACHINE="HPC",OMP_NUM_THREADS=$SIM_THREADS \
        "$HOME/bin/zrun.slurm"
}

submit_snellius() {
    local sim_folder="$1"
    local sim_name="$2"
    
    echo "[$(timestamp)] Submitting (SLURM/Snellius): $sim_name"
    sbatch \
        --job-name="$sim_name" \
        --partition=rome \
        --constraint=scratch-node \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task=$SIM_THREADS \
        --time=5-00:00:00 \
        --output="$sim_folder/%x_%j.out" \
        --error="$sim_folder/%x_%j.err" \
        --export=ALL,SIM_FOLDER="$sim_folder",MACHINE="SNELLIUS",OMP_NUM_THREADS=$SIM_THREADS \
        "$HOME/bin/zrun.slurm"
}

find_simulation_folder() {
    local base_search="$1"
    local sim_name="$2"
    
    # Check base folder
    if [ -d "$base_search" ]; then
        echo "$base_search"
        return 0
    fi
    
    # Check versioned folders
    local i=2
    while [ $i -lt 20 ]; do
        if [ -d "${base_search}_$i" ]; then
            echo "${base_search}_$i"
            return 0
        fi
        i=$((i + 1))
    done
    
    return 1
}


# =====================================================================
# PARSE ARGUMENTS FROM COMMAND LINE
# =====================================================================
NUM_MOLS=""
BASE_DIR=""
MOL_PATH=""

while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help) show_help ;;
        -N|--num_mols) NUM_MOLS="$2"; shift 2 ;;
        -H|--hydronium) H_FILTER="$2"; shift 2 ;;
        -T|--temperatures) T_FILTER="$2"; shift 2 ;;
        -S|--n_steps) N_STEPS="$2"; shift 2 ;;
        -f|--sampling) SAMPLING_FREQ="$2"; shift 2 ;;
        -j|--threads) SIM_THREADS="$2"; shift 2 ;;
        -J|--max_threads) MAX_THREADS="$2"; shift 2 ;;
        -m|--machine) MACHINE="$2"; shift 2 ;;
        --START) RUN_MODE="run"; shift ;;
        --base_dir) BASE_DIR="$2"; shift 2 ;;
        --mol_path) MOL_PATH="$2"; shift 2 ;;
        *)
            echo "[$(timestamp)] Unknown option: $1"
            exit 1
            ;;
    esac
done


# =====================================================================
# VALIDATION
# =====================================================================
if [ -z "$NUM_MOLS" ]; then
    echo "[$(timestamp)] Error: -N/--num_mols is required."
    show_help
fi

MACHINE=$(echo "$MACHINE" | tr '[:lower:]' '[:upper:]')
if [[ ! "$MACHINE" =~ ^(MOPS|HPC|SNELLIUS)$ ]]; then
    echo "[$(timestamp)] Error: Invalid machine type: $MACHINE"
    exit 1
fi


# =====================================================================
# PREPARE PYTHON ARGUMENTS
# =====================================================================
PYTHON_BASE_ARGS="-N $NUM_MOLS -S $N_STEPS -f $SAMPLING_FREQ -j $SIM_THREADS --machine $MACHINE"
if [ -n "$BASE_DIR" ]; then
    PYTHON_BASE_ARGS="$PYTHON_BASE_ARGS --base_dir $BASE_DIR"
    BASE_PATH="$BASE_DIR"
else
    BASE_PATH="$PROJECT_DIR/AMS_simulations"
fi
if [ -n "$MOL_PATH" ]; then
    PYTHON_BASE_ARGS="$PYTHON_BASE_ARGS --mol_path $MOL_PATH"
fi


# =====================================================================
# CONVERT COMMA-SEPARATED LISTS TO ARRAYS
# =====================================================================
IFS=',' read -r -a T_ARR <<< "$T_FILTER"
IFS=',' read -r -a H_ARR <<< "$H_FILTER"


# =====================================================================
# MAIN EXECUTION LOOP
# =====================================================================
echo "[$(timestamp)] Starting simulation setup"
echo "[$(timestamp)] Machine: $MACHINE"
echo "[$(timestamp)] Mode: $RUN_MODE"
echo "[$(timestamp)] Temperatures: [${T_FILTER}]"
echo "[$(timestamp)] H3O counts: [${H_FILTER}]"
echo "[$(timestamp)] Total molecules: $NUM_MOLS"
echo "[$(timestamp)] Threads per sim: $SIM_THREADS"

RUNNING_THREADS=0
TOTAL_SUBMITTED=0
TOTAL_CREATED=0

for T in "${T_ARR[@]}"; do
    for H in "${H_ARR[@]}"; do
        N_H2O=$((NUM_MOLS - H))
        
        if [ $N_H2O -lt 0 ]; then
            echo "[$(timestamp)] Skipping: H3O=$H exceeds total molecules"
            continue
        fi
        
        SIM_NAME="MD_${N_H2O}H2O_${H}H3O_${T}K_${SAMPLING_FREQ}fs"
        LOGFILE="$LOGDIR/${SIM_NAME}.log"
        
        # Create simulation folder
        echo "[$(timestamp)] Creating: $SIM_NAME"
        amspython "$PY" $PYTHON_BASE_ARGS -H $H -T $T > "$LOGFILE" 2>&1
        
        if [ $? -ne 0 ]; then
            echo "[$(timestamp)] FAILED to create: $SIM_NAME"
            continue
        fi
        
        TOTAL_CREATED=$((TOTAL_CREATED + 1))
        
        # Submit if in run mode
        if [ "$RUN_MODE" = "run" ]; then
            # Check thread limit for MOPS/HPC
            if [[ "$MACHINE" =~ ^(MOPS|HPC)$ ]]; then
                if [ $((RUNNING_THREADS + SIM_THREADS)) -gt $MAX_THREADS ]; then
                    echo "[$(timestamp)] Thread limit reached ($RUNNING_THREADS/$MAX_THREADS). Stopping submissions."
                    break 2
                fi
            fi
            
            # Find the created simulation folder
            SIM_FOLDER_BASE="$BASE_PATH/MD_${N_H2O}H2O_${H}H3O/$SIM_NAME"
            SIM_FOLDER=$(find_simulation_folder "$SIM_FOLDER_BASE" "$SIM_NAME")
            
            if [ -z "$SIM_FOLDER" ] || [ ! -d "$SIM_FOLDER" ]; then
                echo "[$(timestamp)] Warning: Could not find folder: $SIM_FOLDER_BASE"
                continue
            fi
            
            # Submit based on machine type
            case "$MACHINE" in
                MOPS)
                    submit_mops "$SIM_FOLDER" "$SIM_NAME"
                    RUNNING_THREADS=$((RUNNING_THREADS + SIM_THREADS))
                    ;;
                HPC)
                    submit_hpc "$SIM_FOLDER" "$SIM_NAME"
                    RUNNING_THREADS=$((RUNNING_THREADS + SIM_THREADS))
                    ;;
                SNELLIUS)
                    submit_snellius "$SIM_FOLDER" "$SIM_NAME"
                    ;;
            esac
            
            TOTAL_SUBMITTED=$((TOTAL_SUBMITTED + 1))
        fi
    done
done


# =====================================================================
# SUMMARY
# =====================================================================
echo ""
echo "[$(timestamp)] ============================================"
echo "[$(timestamp)] Simulation setup completed"
echo "[$(timestamp)] Total folders created: $TOTAL_CREATED"
if [ "$RUN_MODE" = "run" ]; then
    echo "[$(timestamp)] Total simulations submitted: $TOTAL_SUBMITTED"
    if [[ "$MACHINE" =~ ^(MOPS|HPC)$ ]]; then
        echo "[$(timestamp)] Running threads: $RUNNING_THREADS/$MAX_THREADS"
    fi
else
    echo "[$(timestamp)] Test mode - no simulations submitted"
    echo "[$(timestamp)] Use --START to submit simulations"
fi
echo "[$(timestamp)] ============================================"