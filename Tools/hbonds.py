from textwrap import dedent
import matplotlib.pyplot as plt

from base import AnalysisInput, plot_error_handler, log
from common import hist2df, donor_acceptor_cn2df, rdf2df


# =====================================================================
# HBONDS CLASSES
# =====================================================================
class HBonds(AnalysisInput):
    """Parent class for H-bond analysis."""
    FLAG = None
    PLOT_TITLE_PREFIX = "Hydrogen Bonds"
    PLOT_FILENAME_PREFIX = "hbond"


class HBondsWater(HBonds):
    """H-bond analysis for water-water interactions."""
    FLAG = "h_bonds_water"
    CLI_HELP = "Water-water hydrogen bond analysis"
    PLOT_TITLE_PREFIX = "Water-Water Hydrogen Bonds"
    PLOT_FILENAME_PREFIX = "hbond_water"
    
    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        super().__init__(data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps)
        self.filename = ["hist_donor_o_cn", "hist_acceptor_o_cn", "donor_acceptor_cn"]
        self.CONFIG = dedent("""
            HydrogenBond SolventO_H SolventO_H SolventH_O \\
                CopyGroup1 AcceptorO(AllCN) NewGroup1CoordinationGroup 3 \\
                CopyGroup2 DonorO(AllCN) NewGroup2CoordinationGroup 1 \\
                CopyGroup3 DonorH(AllCN) NewGroup3CoordinationGroup 1 

            Histogram {filename_0} DynamicRange GROUPS DonorO(AllCN) PROPERTIES coordinationnumber
            Histogram {filename_1} DynamicRange GROUPS AcceptorO(AllCN) PROPERTIES coordinationnumber
            PrintProperties {filename_2} GROUPS DonorO(AllCN) AcceptorO(AllCN) PROPERTIES timestepiteration coordinationnumber
            """)

    def plot(self):
        """Plot hydrogen bond coordination numbers."""
        self._plot_histogram()
        self._plot_timeseries()
    
    @plot_error_handler("Failed to create H-bond histogram")
    def _plot_histogram(self):
        """Plot histogram of coordination numbers."""
        df_donor = hist2df(self._get_data_file(self.filename[0]))
        df_acceptor = hist2df(self._get_data_file(self.filename[1]))

        plt.figure(figsize=(8, 6))
        plt.bar(df_donor["#value"] - 0.15, df_donor["#fractional"], width=0.3, label="Donors")
        plt.bar(df_acceptor["#value"] + 0.15, df_acceptor["#fractional"], width=0.3, label="Acceptors")
        plt.xlabel("Coordination number")
        plt.ylabel("Fraction")
        plt.title(f"{self.PLOT_TITLE_PREFIX} - {self.sim_tag}")
        plt.legend()
        plt.tight_layout()
        self._save_plot(self.PLOT_FILENAME_PREFIX)
    
    @plot_error_handler("Failed to create H-bond time series")
    def _plot_timeseries(self):
        df = donor_acceptor_cn2df(self._get_data_file(self.filename[2]))
        time_ps = df['timestep'] * self.frame_interval_ps
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        ax1.plot(time_ps, df['avg_donor_cn'], color='blue', alpha=0.7, linewidth=1)
        ax1.set_ylabel('Avg Donor CN')
        ax1.set_title(f'{self.PLOT_TITLE_PREFIX} Time Evolution - {self.sim_tag}')
        ax1.grid(True, alpha=0.3)
        
        # Plot average acceptor coordination
        ax2.plot(time_ps, df['avg_acceptor_cn'], color='red', alpha=0.7, linewidth=1)
        ax2.set_xlabel('Time (ps)')
        ax2.set_ylabel('Avg Acceptor CN')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        self._save_plot(self.PLOT_FILENAME_PREFIX, suffix="_timeseries")


class HBondsMOF(HBonds):
    """H-bond analysis between MOF and water."""
    FLAG = "h_bonds_mof"
    CLI_HELP = "MOF-water hydrogen bond analysis"
    PLOT_TITLE_PREFIX = "MOF-Water Hydrogen Bonds"
    PLOT_FILENAME_PREFIX = "hbond_mof"
    
    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        super().__init__(data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps)
        self.filename = ["hbond_mu3o_h", "hbond_mu3h_o"]
        self.CONFIG = dedent("""
            # MOF-Solvent H-bonds
            HydrogenBond Mu3O_H SolventO_H SolventH_O Filename {filename_0}
            HydrogenBond SolventO_H Mu3O_H Mu3H_O Filename {filename_1}
            """)
    
    def get_output_patterns(self):
        return [f"{self.output_dir}{fname}_4_{{tag}}.txt" for fname in self.filename]
    
    @plot_error_handler("Failed to create MOF H-bond RDF")
    def plot(self):
        """Plot MOF-water hydrogen bond RDFs."""
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        labels = ["μ3-O···Water", "Linker-O···Water"]
        colors = ['red', 'blue']
        
        for idx, (fname, label, color) in enumerate(zip(self.filename, labels, colors)):
            # Read O...O RDF (file suffix _4 from HydrogenBond output)
            rdf_file = self._get_data_file(f"{fname}_4")
            if rdf_file.exists():
                df = rdf2df(rdf_file)
                bond_col = df.columns[0]  # distance
                weight_col = df.columns[1]  # normalized weight
                axes[idx].plot(df[bond_col], df[weight_col], color=color, linewidth=1.5)
                axes[idx].set_xlabel("O···O Distance (Å)")
                axes[idx].set_ylabel("g(r)")
            else:
                axes[idx].text(0.5, 0.5, f"No data: {rdf_file.name}", 
                              ha='center', va='center', transform=axes[idx].transAxes)
        
        plt.suptitle(f"{self.PLOT_TITLE_PREFIX} - {self.sim_tag}")
        plt.tight_layout()
        self._save_plot(self.PLOT_FILENAME_PREFIX)
