from textwrap import dedent
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from base import AnalysisInput, plot_error_handler, log
from common import (
    extract_ptfel_barrier,
)


# =====================================================================
# PTFEL CLASSES
# =====================================================================
class Ptfel(AnalysisInput):
    """Parent class for PTFEL analysis."""
    FLAG = None
    PLOT_TITLE_PREFIX = "Proton Transfer Free Energy Landscape"
    PLOT_FILENAME_PREFIX = "ptfel"
    XLABEL_1D = "δ = d(O*-H) - d(H-O) (Å)"
    YLABEL_2D = "O*-O Distance (Å)"
    CMAP_2D = "viridis_r"
    LINE_COLOR = 'b'
    
    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        super().__init__(data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps)
        self.delta_min = np.nan
        self.delta_saddle = np.nan
        self.delta_1d_limit = 2.0  # Å
        self.spline_func = None  # Smoothing spline or polynomial fit function
    
    def post_process(self):
        """Rename 2D output file."""
        old_path = self.data_dir / f"{self.filename[0]}_{self.sim_tag}.txt2d"
        new_path = self.data_dir / f"{self.filename[1]}_{self.sim_tag}.txt"
        try:
            if old_path.exists():
                if new_path.exists():
                    new_path.unlink()  # Delete old destination file first
                old_path.rename(new_path)  # Then rename new source to destination
                log(f"[PostProcess] Renamed: {old_path.name} -> {new_path.name}")
        except Exception as e: 
            log(f"[PostProcess] Error: {e}")
    
    def plot(self):
        """Plot proton transfer free energy landscapes (1D and 2D)."""
        self._plot_1d()
        self._plot_2d()
    
    @plot_error_handler()
    def _plot_1d(self):
        """Plot 1D free energy landscape with barrier extraction."""
        ptfel_file = self._get_data_file(self.filename[0])
        df = pd.read_csv(ptfel_file, sep=r"\s+", header=None, names=["delta", "neg_ln_count", "count"])
        df_valid = df[df["neg_ln_count"] != np.inf].copy()
        
        if df_valid.empty:
            log(f"[Plotter] {self.PLOT_FILENAME_PREFIX} 1D: No proton transfer events detected")
            return
        
        delta = df_valid["delta"].to_numpy()
        fe_raw = df_valid["neg_ln_count"].to_numpy()
        fe_shifted = fe_raw - fe_raw.min()
        
        # Extract barrier using SavGol smoothing
        result = extract_ptfel_barrier(delta, fe_shifted)
        self.delta_saddle = result['delta_saddle']
        self.fe_saddle = result['fe_saddle']
        self.delta_min = result['delta_min']
        self.spline_func = result['spline_func']
        
        # Store raw data and shifted data for comparison plots
        self.delta_raw = delta
        self.fe_raw = fe_shifted
        
        # Shift delta so saddle is at origin
        delta_centered = delta - self.delta_saddle if not np.isnan(self.delta_saddle) else delta
        self.delta_centered = delta_centered
        
        # Plot with centered delta (saddle at origin)
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(delta_centered, fe_shifted, f'{self.LINE_COLOR}-', linewidth=1, alpha=0.4, label='Raw data')
        
        # Plot SavGol smoothed fit if available
        if self.spline_func is not None:
            delta_fine = np.linspace(delta.min(), delta.max(), 300)
            delta_fine_centered = delta_fine - self.delta_saddle if not np.isnan(self.delta_saddle) else delta_fine
            ax.plot(delta_fine_centered, self.spline_func(delta_fine), 'k-', linewidth=2.5, label='SavGol smoothed')
            
        # Saddle is now at origin (delta=0)
        ax.axvline(0, color='r', linestyle='--', linewidth=1.5, alpha=0.7, label='Saddle (delta=0)')
        # Mark saddle point (barrier) at origin
        ax.scatter([0], [self.fe_saddle], color='r', s=120, zorder=5, marker='^', edgecolors='k', linewidth=1)
        
        # Paper-quality styling
        ax.set_xlabel("delta (Å)", fontsize=14, fontweight='bold')
        ax.set_ylabel("Free Energy (kT)", fontsize=14, fontweight='bold')
        ax.set_xlim(-self.delta_1d_limit, self.delta_1d_limit)
        ax.set_ylim(-0.2, max(fe_shifted.max(), self.fe_saddle + 0.5) + 0.5)
        ax.set_title(f"{self.PLOT_TITLE_PREFIX} - {self.sim_tag}", fontsize=16, fontweight='bold')
        ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
        ax.grid(True, alpha=0.4, linestyle='-', linewidth=0.5)
        ax.tick_params(axis='both', which='major', labelsize=12)
        plt.tight_layout()
        self._save_plot(self.PLOT_FILENAME_PREFIX)
    
    @plot_error_handler()
    def _plot_2d(self):
        """Plot 2D free energy landscape."""
        ptfel2d_file = self._get_data_file(self.filename[1])
        df = pd.read_csv(ptfel2d_file, sep=r"\s+", header=None, names=["delta", "distance", "free_energy", "count"])
        df_valid = df[df["free_energy"] != np.inf].copy()
        
        if df_valid.empty:
            log(f"[Plotter] {self.PLOT_FILENAME_PREFIX} 2D: No valid data points")
            return
        
        df_valid["free_energy"] = df_valid["free_energy"] - df_valid["free_energy"].min()
        pivot = df_valid.pivot_table(index="distance", columns="delta", values="free_energy", aggfunc="mean")
        pivot = pivot.sort_index(ascending=True).reindex(sorted(pivot.columns), axis=1)
        
        plt.figure(figsize=(10, 8))
        delta_vals = pivot.columns.values
        dist_vals = pivot.index.values
        
        delta_edges = np.concatenate([[delta_vals[0] - (delta_vals[1] - delta_vals[0]) / 2],
            (delta_vals[:-1] + delta_vals[1:]) / 2,
            [delta_vals[-1] + (delta_vals[-1] - delta_vals[-2]) / 2]])
        dist_edges = np.concatenate([[dist_vals[0] - (dist_vals[1] - dist_vals[0]) / 2],
            (dist_vals[:-1] + dist_vals[1:]) / 2,
            [dist_vals[-1] + (dist_vals[-1] - dist_vals[-2]) / 2]])
        
        im = plt.pcolormesh(delta_edges, dist_edges, pivot.values, cmap=self.CMAP_2D, shading="flat")
        plt.colorbar(im, label="Free Energy (kT)")
        plt.xlabel(self.XLABEL_1D)
        plt.ylabel(self.YLABEL_2D)
        plt.xlim(-self.delta_1d_limit, self.delta_1d_limit)
        plt.ylim(2.2, 3.5)
        plt.title(f"2D {self.PLOT_TITLE_PREFIX} - {self.sim_tag}")
        plt.tight_layout()
        self._save_plot(self.PLOT_FILENAME_PREFIX, suffix="2d")


class PtfelWater(Ptfel):
    """PTFEL for water-water proton transfer."""
    FLAG = "ptfel_water"
    CLI_HELP = "Water-water proton transfer analysis"
    PLOT_TITLE_PREFIX = "Water-Water PTFEL"
    PLOT_FILENAME_PREFIX = "ptfel_water"
    
    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        super().__init__(data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps)
        self.filename = ["ptfel_water", "ptfel_water_2d"]
        # Water-water PTFEL using coordinated groups
        self.CONFIG = dedent("""
            # Water-water h-bond groups (for water PTFEL)
            HydrogenBond SolventO_H SolventO_H SolventH_O \\
                CopyGroup1 AcceptorO NewGroup1CoordinationGroup 3 \\
                CopyGroup2 DonorO NewGroup2CoordinationGroup 1 \\
                CopyGroup3 DonorH NewGroup3CoordinationGroup 1

            # Water-water proton transfer: AcceptorO (hydronium) ← H → DonorO (water)
            DoubleCoordinationShortDelta Filename {filename_0} \\
                LHS AcceptorO SolventH_O LGroup3MustBe DonorO \\
                RHS AcceptorO SolventH_O RGroup3MustBe DonorO \\
                Every 1 MinValue 0 MaxValue 2.0 Resolution 0.01 \\
                Min13 2.0 Max13 3.5 Resolution13 0.01 WellAtZero
            """)


class PtfelMOF(Ptfel):
    """PTFEL for MOF-water proton transfer."""
    FLAG = "ptfel_mof"
    CLI_HELP = "MOF-water proton transfer analysis"
    PLOT_TITLE_PREFIX = "MOF-Water PTFEL"
    PLOT_FILENAME_PREFIX = "ptfel_mof"
    XLABEL_1D = "delta = d(MOF-O-H) - d(H-O_water) (Å)"
    YLABEL_2D = "MOF-O···O Distance (Å)"
    CMAP_2D = "magma_r"
    LINE_COLOR = 'r'
    
    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        super().__init__(data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps)
        self.filename = ["ptfel_mof", "ptfel_mof_2d"]
        # MOF-water PTFEL using asymmetric DCSD
        # LHS (delta>0): H closer to water donor, MOF is acceptor
        # RHS (delta<0): H closer to MOF donor, water is acceptor
        self.CONFIG = dedent("""
            # MOF-Solvent h-bond groups (for MOF PTFEL)

            # H-bond: MOF O (acceptor) < H < Solvent O (donor)
            HydrogenBond MofO_H SolventO_H SolventH_O \\
                CopyGroup1 MofAcceptorO NewGroup1CoordinationGroup 3 \\
                CopyGroup2 SolventDonorO NewGroup2CoordinationGroup 1 \\
                CopyGroup3 SolventDonorH NewGroup3CoordinationGroup 1

            # H-bond: Solvent O (acceptor) < H < MOF O (donor)
            HydrogenBond SolventO_H MofO_H MofH_O \\
                CopyGroup1 SolventAcceptorO NewGroup1CoordinationGroup 3 \\
                CopyGroup2 MofDonorO NewGroup2CoordinationGroup 1 \\
                CopyGroup3 MofDonorH NewGroup3CoordinationGroup 1

            # Create combined H group with O coordination for DCSD Group2
            DefineGroup AllH SUM SolventH MofH
            DefineGroup AllO SUM MofO SolventO
            DefineGroup AllO_H FINDSHORTEST FromGroup AllH ToGroup AllO
            DefineGroup AllH_O INVERTCOORDINATION AllO_H

            # CRITICAL: For DCSD, Group1 needs coordination matching Group2 (AllH_O)
            # Use INTERSECTION+SUM pattern to give acceptor groups AllO_H coordination
            DefineGroup MofAcceptorO_AllH INTERSECTION AllO_H MofAcceptorO
            ModifyGroup MofAcceptorO_AllH SUM MofAcceptorO_AllH MofAcceptorO

            DefineGroup SolventAcceptorO_AllH INTERSECTION AllO_H SolventAcceptorO
            ModifyGroup SolventAcceptorO_AllH SUM SolventAcceptorO_AllH SolventAcceptorO

            # MOF-Solvent proton transfer: MOF_O < H > SolventO
            # LHS (delta>0): H closer to water donor, MOF is acceptor
            # RHS (delta<0): H closer to MOF donor, water is acceptor
            DoubleCoordinationShortDelta Filename {filename_0} \\
                LHS MofAcceptorO_AllH AllH_O LGroup3MustBe SolventDonorO \\
                RHS SolventAcceptorO_AllH AllH_O RGroup3MustBe MofDonorO \\
                Every 1 MinValue 0.0 MaxValue 2 Resolution 0.01 \\
                Min13 2.0 Max13 3.5 Resolution13 0.01 WellAtZero
            """)
