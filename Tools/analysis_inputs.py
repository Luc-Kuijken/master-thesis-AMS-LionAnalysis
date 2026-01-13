#!/usr/bin/env amspython
# =====================================================================
# analysis_inputs.py
# ---------------------------------------------------------------------
# Analysis input definitions for LionAnalysis trajectory analysis.
# 
# This module re-exports:
#   - AnalysisInput base class
#   - log() and plot_error_handler() utilities
#   - AnalysisRegistry for managing analysis instances
#   - All child analysis classes
#
# Child classes:
#   - rdf.py: Rdf, RdfWater, RdfMOF
#   - hbonds.py: HBonds, HBondsWater, HBondsMOF
#   - ptfel.py: Ptfel, PtfelWater, PtfelMOF
#   - msd.py: Msd, MsdWater, MsdHydronium
#   - residence_times.py: ResidenceTimes
#   - debug.py: Debug, AtomXyz, TrajectorySnapshot
#
# Author: L. Kuijken
# Last updated: 2025-12-13
# =====================================================================

# Import base class and utilities from base.py
from base import AnalysisInput, log, plot_error_handler


# =====================================================================
# REGISTRY
# =====================================================================
class AnalysisRegistry: 
    """Central registry for managing analysis instances."""
    
    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        self.data_dir = data_dir
        self.msd_dir = msd_dir
        self.plots_dir = plots_dir
        self.sim_tag = sim_tag
        self.frame_interval_ps = frame_interval_ps
        self._instances = {}
        self._initialize_all()
    
    def _initialize_all(self):
        for flag, cls in AnalysisInput.get_all_analysis_classes().items():
            self._instances[flag] = cls(
                self.data_dir, self.msd_dir, self.plots_dir, 
                self.sim_tag, self.frame_interval_ps
            )
    
    def get(self, flag):
        return self._instances.get(flag)
    
    def get_all(self):
        return self._instances
    
    def __getitem__(self, flag):
        return self.get(flag)
    
    def items(self):
        return self._instances.items()


# =====================================================================
# IMPORT CHILD CLASSES
# =====================================================================
# RDF analysis
from rdf import Rdf, RdfWater, RdfMOF
# Hydrogen bond analysis
from hbonds import HBonds, HBondsWater, HBondsMOF
# Proton transfer free energy landscape
from ptfel import Ptfel, PtfelWater, PtfelMOF
# Mean square displacement
from msd import Msd, MsdWater, MsdHydronium
# Residence times
from residence_times import ResidenceTimes
# Debug and utility classes
from debug import Debug, AtomXyz, TrajectorySnapshot
