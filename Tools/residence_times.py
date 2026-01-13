from textwrap import dedent
import numpy as np
import matplotlib.pyplot as plt

from base import AnalysisInput, plot_error_handler, log
from common import (
    tcf2df, extract_lifetime_from_tcf,
)


# =====================================================================
# RESIDENCE TIMES
# =====================================================================
class ResidenceTimes(AnalysisInput):
    """Residence time analysis using SSP."""
    FLAG = "res_times"
    CLI_HELP = "Residence time analysis (SSP)"
    FIT_METHOD = "double_exp"
    
    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        super().__init__(data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps)
        self.filename = ["resid_hydronium"]
        self.lifetimes = {}
        self.hopping_rates = {}
        self.CONFIG = dedent("""
            TResidenceTimeSSP ParentGroup SolventO_H \\
                Reactants HydroniumO_H Products WaterO_H Filename {filename_0} \\
                OldValueEscape1 0 0 0 UR NewValueEscape2! 0 0 0 VR \\
                MaxHistory {max_history} RealTime PrintEvery 10
            """)
    
    def plot(self):
        for fname in self.filename:
            self._plot_single_residence(fname)
        if self.lifetimes:
            log(f"[ResidenceTimes] Extracted lifetimes:")
            for name, tau in self.lifetimes.items():
                log(f"[ResidenceTimes]   {name}: tau = {tau:.3f} ps")
    
    @plot_error_handler("Failed to create residence plot")
    def _plot_single_residence(self, fname):
        f = self._get_data_file(fname)
        df = tcf2df(f)
        species = "hydronium" if "hydronium" in fname else "water"
        
        if df.empty:
            log(f"[ResidenceTimes] No data in {fname}")
            return
        
        t = df["#t(ps)"].to_numpy()
        y = df["value"].to_numpy()
        y_norm = y / y[0] if y[0] != 0 else y
        
        tau, fit_info = extract_lifetime_from_tcf(t, y_norm, method=self.FIT_METHOD)
        
        # Store extracted values
        self.lifetimes[species] = tau
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(t, y_norm, 'b-', linewidth=1, label="C(t) data", alpha=0.7)
        
        if "fitted" in fit_info:
            ax.plot(t, fit_info["fitted"], 'r--', linewidth=2, label=f"Fit ({self.FIT_METHOD})")
        
        if not np.isnan(tau):
            textstr = f"tau = {tau:.3f} ps"
            if self.FIT_METHOD == "double_exp" and "A" in fit_info:
                textstr += f"\nA = {fit_info['A']:.3f}, tau1 = {fit_info['tau1']:.3f}, tau2 = {fit_info['tau2']:.3f}"
            ax.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                   verticalalignment='top', horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel("C(t)")
        ax.set_title(f"Residence Time - {species.capitalize()} ({self.sim_tag})")
        ax.legend()
        ax.set_xlim(0, t.max())
        ax.set_ylim(0, 1.1)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        self._save_plot(f"resid_{species}")

