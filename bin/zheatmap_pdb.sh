#!/bin/bash
# =====================================================================
# zheatmap_pdb.sh
# ---------------------------------------------------------------------
# Parallel heatmap analysis launcher for SITES-ANALYZER.
# Processes PDB trajectories to generate Average Occupation Profiles.
#
# Usage:
#   zheatmap_pdb.sh [BASE_FOLDER]
#
# Author: L. Kuijken
# Last updated: 2025-11-19
# =====================================================================


# =====================================================================
# INITIAL SETUP
# =====================================================================
TRAJ_DIR="${1:-$HOME/master_thesis_project/Trajectories/pdb/}"
ANALYSIS_DIR="$HOME/master_thesis_project/SitesAnalysis"
SITES_CONFIG="$HOME/bin/sites.input"
LOGDIR="$ANALYSIS_DIR/Analysis_logs"
MASTERLOG="$LOGDIR/master_summary.log"
DATETIME=$(date '+%Y-%m-%d %H:%M:%S')
TIME=$(date +%H:%M:%S)

mkdir -p "$LOGDIR"
echo "" >> "$MASTERLOG"
echo "[$DATETIME]: Running zheatmap_pdb.sh > $MASTERLOG..."
find "$LOGDIR" -type f ! -name '*master*.log' -delete

# =====================================================================
# BACKGROUND EXECUTION
# =====================================================================
# Execute in background with all output to logs
if [ -z "$HEATMAP_BACKGROUND" ]; then
    export HEATMAP_BACKGROUND=1
    exec "$0" "$@" >> "$MASTERLOG" 2>&1 &
    exit $?
fi


# =====================================================================
# VALIDATION
# =====================================================================
if [ ! -d "$TRAJ_DIR" ]; then
    echo "[$TIME] Error: Trajectory directory not found: $TRAJ_DIR"
    exit 1
fi


# =====================================================================
# TRAJECTORY DISCOVERY
# =====================================================================
JOBLOG="$LOGDIR/heatmap_$(date +%Y-%m-%d_%H:%M:%S).log"

readarray -t FILES < <(find "$TRAJ_DIR" -type f -name "*_traj.pdb" 2>/dev/null | sort -u)
TOTAL=${#FILES[@]}

if [ "$TOTAL" -eq 0 ]; then
    echo "[$TIME] No .pdb trajectories found in $TRAJ_DIR" | tee -a "$JOBLOG"
    exit 0
fi


# =====================================================================
# PARALLEL EXECUTION 
# =====================================================================
echo "[$TIME] Found $TOTAL PDB trajectories in $TRAJ_DIR" | tee -a "$JOBLOG"
echo "[$TIME] Command: $0 $TRAJ_DIR"
echo "[$TIME] Heatmap analysis started"

COUNT=0
FAILED=0

for pdb in "${FILES[@]}"; do
    COUNT=$((COUNT + 1))
    basename=$(basename "$pdb" _traj.pdb)
    LOGFILE="$LOGDIR/${basename}_heatmap.log"
    
    echo "[$TIME] [$COUNT/$TOTAL] Launching: $basename" | tee -a "$JOBLOG"
    
    # Launch analysis in background
    (
        bash ~/bin/zzrunSITES "$pdb" "$SITES_CONFIG" "$ANALYSIS_DIR" >> "$LOGFILE" 2>&1
        STATUS=$?
        if [ $STATUS -eq 0 ]; then
            echo "[$TIME] [$COUNT/$TOTAL] Completed: $basename" | tee -a "$JOBLOG"
        else
            echo "[$TIME] [$COUNT/$TOTAL] FAILED: $basename (exit $STATUS)" | tee -a "$JOBLOG"
            echo "$basename" >> "$LOGDIR/failed_heatmaps.txt"
        fi
    ) &
done

# Wait for all background jobs to complete
wait

# =====================================================================
# SUMMARY
# =====================================================================
if [ -f "$LOGDIR/failed_heatmaps.txt" ]; then
    FAILED=$(wc -l < "$LOGDIR/failed_heatmaps.txt")
fi

SUCCESSFUL=$((TOTAL - FAILED))

echo "===================================================================="
echo "[$TIME] Heatmap analysis completed"
echo "[$TIME] Total trajectories: $TOTAL"
echo "[$TIME] Successful: $SUCCESSFUL"
echo "[$TIME] Failed: $FAILED"
echo "===================================================================="

if [ "$FAILED" -gt 0 ]; then
    echo "[$TIME] Failed analyses listed in: $LOGDIR/failed_heatmaps.txt"
fi