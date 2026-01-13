from textwrap import dedent
import matplotlib.pyplot as plt

from base import AnalysisInput, plot_error_handler, log
from common import rdf2df


# =====================================================================
# RDF CLASSES
# =====================================================================
class Rdf(AnalysisInput):
    """Parent class for RDF analysis."""
    FLAG = None
    PLOT_TITLE_PREFIX = "Radial Distribution Function"
    PLOT_FILENAME_PREFIX = "rdf"
    
    @plot_error_handler("Failed to create RDF plot")
    def plot(self):
        plt.figure(figsize=(8, 6))
        for fname in self.filename:
            f = self._get_data_file(fname)
            df = rdf2df(f)
            if "#r" in df.columns and "RDF" in df.columns:
                label = fname.replace("rdf_", "").replace("_", "-")
                plt.plot(df["#r"], df["RDF"], label=label)
        plt.xlabel("r (Å)")
        plt.ylabel("g(r)")
        plt.title(f"{self.PLOT_TITLE_PREFIX} - {self.sim_tag}")
        plt.legend()
        plt.tight_layout()
        self._save_plot(self.PLOT_FILENAME_PREFIX)


class RdfWater(Rdf):
    """RDF for H2O-H2O interactions."""
    FLAG = "rdf_water"
    CLI_HELP = "RDF for H2O-H2O interactions"
    PLOT_TITLE_PREFIX = "H2O-H2O Radial Distribution Function"
    PLOT_FILENAME_PREFIX = "rdf_water"
    
    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        super().__init__(data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps)
        self.filename = ["rdf_o_h", "rdf_o_o"]
        self.CONFIG = dedent("""
            # --- RDFs Solvent-Solvent ---
            RDF SolventO_H SolventH_O {filename_0} MinDist 0.0 MaxDist 8.0 Resolution 0.01
            RDF SolventO_H SolventO_H {filename_1} MinDist 0.0 MaxDist 8.0 Resolution 0.01
            """)


class RdfMOF(Rdf):
    """RDF for MOF-H2O interactions."""
    FLAG = "rdf_mof"
    CLI_HELP = "RDF for MOF-H2O interactions"
    PLOT_TITLE_PREFIX = "H2O-MOF Radial Distribution Function"
    PLOT_FILENAME_PREFIX = "rdf_mof"
    
    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        super().__init__(data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps)
        self.filename = ["rdf_mu3o_h", "rdf_mu3o_o", "rdf_linkero_h", "rdf_linkero_o"]
        self.CONFIG = dedent("""
            # --- RDFs Mu3O-Solvent ---
            RDF Mu3O_H SuperH_O {filename_0} MinDist 0.0 MaxDist 8.0 Resolution 0.01
            RDF Mu3O_H SolventO_H {filename_1} MinDist 0.0 MaxDist 8.0 Resolution 0.01
                    
            # --- RDFs LinkerO-Solvent ---
            RDF LinkerO_CH SuperH_O {filename_2} MinDist 0.0 MaxDist 8.0 Resolution 0.01
            RDF LinkerO_CH SolventO_H {filename_3} MinDist 0.0 MaxDist 8.0 Resolution 0.01
            """)