#!/usr/bin/env amspython
# =====================================================================
# base.py
# ---------------------------------------------------------------------
# Base class for LionAnalysis trajectory analysis.
#
# Author: L. Kuijken
# Last updated: 2025-12-13
# =====================================================================

from pathlib import Path
from datetime import datetime
from textwrap import dedent
from functools import wraps
import matplotlib.pyplot as plt


def log(message):
    """Print message with timestamp."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {message}", flush=True)


def plot_error_handler(error_message=None):
    """Decorator for plot functions with error handling."""
    def decorator(plot_func):
        @wraps(plot_func)
        def wrapper(self, *args, **kwargs):
            try:
                return plot_func(self, *args, **kwargs)
            except Exception as e:
                msg = error_message or f"Failed to create {plot_func.__name__}"
                log(f"[Plotter] {msg}: {e}")
        return wrapper
    return decorator


# =====================================================================
# BASE CLASS
# =====================================================================
class AnalysisInput:
    """Base class for all analysis input configurations."""

    # Class registry for all analysis types
    _registry = {}

    # Override these in subclasses
    FLAG = None
    CLI_HELP = "No description"
    CLI_ACTION = "store_true"
    CLI_VAR = "VALUE"

    # Common configuration template (class attribute)
    COMMON_CONFIG = dedent(
        """
        DumpFile {xyz_file}
        DumpFileFormat mattixyz
        Threads {n_threads}
        Overwrite
        CoutFrequency 500
        IntelligentUnwrap
        BasicTimeUnit {timestep_ps}
        SuperEvery {stride_length}
        Prefix {output_dir}/
        Suffix _{sim_tag}.txt
        """)

    GROUP_DEF_CONFIG = dedent(
        """
        # --- Atom Groups ---
        DefineGroup H ATOMICNUMBER H
        DefineGroup O ATOMICNUMBER O
        DefineGroup Zr ATOMICNUMBER Zr
        DefineGroup C ATOMICNUMBER C

        # --- Solvent and MOF Groups ---
        DefineGroup Solvent LIST {solvent_indices}
        DefineGroup Mof DIFF All Solvent

        # --- Coordinated O Groups ---
        # SuperO: O atoms with their closest H (global coordination)
        DefineGroup SuperO FINDSHORTEST FromGroup H ToGroup O
        DefineGroup O_H BOND O H MaxDist 1.25
        ModifyGroup SuperO INTERSECTION SuperO O_H

        # Solvent O with H coordination
        DefineGroup SolventO INTERSECTION O Solvent
        DefineGroup SolventO_H INTERSECTION SuperO SolventO
        ModifyGroup SolventO_H SUM SolventO_H SolventO
        
        # MOF O with H coordination
        DefineGroup MofO INTERSECTION O Mof
        DefineGroup MofO_H INTERSECTION SuperO MofO
        ModifyGroup MofO_H SUM MofO_H MofO

        # --- Coordinated H Groups ---
        DefineGroup SolventH INTERSECTION H Solvent
        DefineGroup SuperH_O INVERTCOORDINATION SuperO
        DefineGroup SolventH_O INVERTCOORDINATION SolventO_H
        DefineGroup MofH INTERSECTION H Mof
        DefineGroup MofH_O INVERTCOORDINATION MofO_H

        # --- Coordinated Solvent O Groups (by CN) ---
        DefineGroup HydroniumO_H SUBGROUPCN SolventO_H Coord 3
        DefineGroup WaterO_H SUBGROUPCN SolventO_H Coord 2
        DefineGroup HydroxideO_H SUBGROUPCN SolventO_H Coord 1
        DefineGroup RadicalO SUBGROUPCN SolventO_H Coord 0

        # --- Coordinated Solvent H Groups ---
        DefineGroup HydroniumH_O INVERTCOORDINATION HydroniumO_H
        DefineGroup WaterH_O INVERTCOORDINATION WaterO_H
        DefineGroup HydroxideH_O INVERTCOORDINATION HydroxideO_H
        DefineGroup RadicalH SUBGROUPCN SolventH_O Coord 0

        # --- Coordinated linker Groups ---
        DefineGroup LinkerO_CH BOND MofO_H C MaxDist 1.6
        DefineGroup LinkerCH_O INVERTCOORDINATION LinkerO_CH
        DefineGroup LinkerC_O INTERSECTION LinkerCH_O C
        DefineGroup LinkerH_O INTERSECTION LinkerCH_O H

        # --- Coordinated mu3 Groups ---
        DefineGroup Mu3O_H DIFF MofO_H LinkerO_CH
        DefineGroup Mu3H_O INVERTCOORDINATION Mu3O_H
        
        DefineGroup TrackedAtom LIST {atom_indices} Static
        """)

    def __init_subclass__(cls, **kwargs):
        """Automatically register all subclasses."""
        super().__init_subclass__(**kwargs)
        if cls.FLAG is not None:
            AnalysisInput._registry[cls.FLAG] = cls

    @classmethod
    def get_all_analysis_classes(cls):
        """Return all registered analysis classes."""
        return cls._registry

    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        self.data_dir = Path(data_dir)
        self.msd_dir = Path(msd_dir)
        self.plots_dir = Path(plots_dir)
        self.sim_tag = sim_tag
        self.frame_interval_ps = frame_interval_ps
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self.flag = self.FLAG
        
        # To be set by child classes
        self.CONFIG = ""
        self.filename = []
        self.output_dir = "data/"
    
    def _get_data_file(self, filename_base, output_dir=None):
        """Helper to construct full data file path."""
        if output_dir is None:
            output_dir = self.data_dir
        return Path(output_dir) / f"{filename_base}_{self.sim_tag}.txt"
    
    def _save_plot(self, filename_prefix, suffix="", dpi=150):
        """Helper to save plot with consistent naming."""
        filename = f"{filename_prefix}{suffix}_{self.sim_tag}.png"
        outpath = self.plots_dir / filename
        plt.savefig(outpath, dpi=dpi)
        plt.close()
        log(f"[Plotter] Saved: {outpath.name}")
        return outpath
    
    def get_output_patterns(self):
        """Return list of output file patterns with {tag} placeholder."""
        return [f"{self.output_dir}{fname}_{{tag}}.txt" for fname in self.filename]

    def plot(self):
        """Generate plots for this analysis. Override in child classes."""
        pass
    
    def post_process(self):
        """Post-processing step after analysis. Override in child classes if needed."""
        pass
