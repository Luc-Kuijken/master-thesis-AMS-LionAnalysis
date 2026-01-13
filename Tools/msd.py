from textwrap import dedent
import matplotlib.pyplot as plt

from base import AnalysisInput, plot_error_handler, log
from common import msd2df


# =====================================================================
# MSD CLASSES
# =====================================================================
class Msd(AnalysisInput):
    """Base class for MSD analysis."""
    FLAG = None
    PLOT_TITLE_PREFIX = "Mean Square Displacement"
    PLOT_FILENAME_PREFIX = "msd"
    YLABEL = "MSD (Å²)"
    PLOT_LABEL = "MSD"
    
    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        super().__init__(data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps)
        self.output_dir = "msd/"
    
    @plot_error_handler("Failed to create MSD plot")
    def plot(self):
        msd_file = self._get_data_file(self.filename[0], output_dir=self.msd_dir)
        df = msd2df(msd_file)
        
        plt.figure(figsize=(8, 6))
        plt.plot(df["#t(ps)"], df["value"], label=self.PLOT_LABEL)
        plt.xlabel("Time (ps)")
        plt.ylabel(self.YLABEL)
        plt.title(f"{self.PLOT_TITLE_PREFIX} - {self.sim_tag}")
        plt.legend()
        plt.tight_layout()
        self._save_plot(self.PLOT_FILENAME_PREFIX)


class MsdWater(Msd):
    """MSD for water."""
    FLAG = "msd_water"
    CLI_HELP = "Water MSD analysis"
    PLOT_TITLE_PREFIX = "Mean Square Displacement (Water)"
    PLOT_FILENAME_PREFIX = "msd_water"
    PLOT_LABEL = "water"
    
    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        super().__init__(data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps)
        self.filename = ["msd_water"]
        self.CONFIG = dedent("""
            TMSD ParentGroup SolventO_H Group SolventO_H Filename {filename_0} MaxHistory {max_history} PrintEvery 1
            """)


class MsdHydronium(Msd):
    """MSD for hydronium (vehicular diffusion)."""
    FLAG = "msd_hydronium"
    CLI_HELP = "Hydronium MSD analysis"
    PLOT_TITLE_PREFIX = "Mean Square Displacement (Hydronium)"
    PLOT_FILENAME_PREFIX = "msd_hydronium"
    PLOT_LABEL = "hydronium O"
    
    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        super().__init__(data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps)
        self.filename = ["msd_hydronium_o"]
        # Vehicular MSD: track H3O+ oxygens with MemberEscape
        # NewMemberEscape1 -1 5 -1 UR @: exclude from statistics if not continuously H3O+ for 5 frames
        self.CONFIG = dedent("""
            TMSDFollow ParentGroup SolventO_H Group HydroniumO_H Filename {filename_0} \\
                MaxHistory {max_history} RealTime PrintEvery 500 \\
                NewMemberEscape1 -1 5 -1 VR @
            """)
