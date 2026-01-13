from textwrap import dedent
import numpy as np
import matplotlib.pyplot as plt

from base import AnalysisInput, plot_error_handler, log
from common import averagesize2df


# =====================================================================
# DEBUG AND OTHER
# =====================================================================
class Debug(AnalysisInput):
    """Debug diagnostics."""
    FLAG = "debug"
    CLI_HELP = "Debug diagnostics"
    
    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        super().__init__(data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps)
        self.filename = ["averagesizes", "groupsizes_framewise"]
        self.CONFIG = dedent("""
            FinalPrintProperties {filename_0} MultipleLines PrintEvery 100 ALLGROUPS PROPERTIES groupname numatoms averagesize totalsize
            PrintProperties {filename_1} Every 1 ALLGROUPS PROPERTIES groupname numatoms
            """)
    
    @plot_error_handler("Failed to create average size plot")
    def plot(self):
        df = averagesize2df(self._get_data_file(self.filename[0]))
        plt.figure(figsize=(8, 6))
        plt.barh(df["Group"], df["AverageSize"])
        plt.xlabel("Average group size")
        plt.ylabel("Group")
        plt.title(f"Group Size Diagnostics - {self.sim_tag}")
        plt.tight_layout()
        self._save_plot("averagesize")


class AtomXyz(AnalysisInput):
    """Atom trajectory tracking."""
    FLAG = "atom_xyz"
    CLI_HELP = "Track atoms by indices"
    CLI_ACTION = "store"
    CLI_VAR = "INDICES"
    
    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        super().__init__(data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps)
        self.filename = ["atom_coords"]
        self.CONFIG = dedent("""
            PrintProperties {filename_0} Every 100 TabSeparator GROUPS TrackedAtom PROPERTIES timestepiteration id type x y z
            """)
    
    @plot_error_handler("Failed to create atom trajectory plot")
    def plot(self):
        atom_coords_file = self._get_data_file(self.filename[0])
        with open(atom_coords_file, 'r') as file:
            lines = [line.strip() for line in file if not line.startswith('#')]
        
        if not lines:
            return
        
        first_line = lines[0].split()
        n_atoms = (len(first_line) - 1) // 5
        if n_atoms == 0:
            return
        
        # Extract atom info
        atom_info = []
        for i in range(n_atoms):
            offset = 1 + i * 5
            atom_info.append((int(first_line[offset]), first_line[offset + 1]))
        
        # Load data
        usecols = [0]
        for i in range(n_atoms):
            id_col = 1 + i * 5
            usecols.extend([id_col, id_col + 2, id_col + 3, id_col + 4])
        
        data = np.loadtxt(atom_coords_file, usecols=usecols)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        
        color_map = {'O': '#FF4444', 'H': '#4444FF', 'C': '#888888', 'N': '#44FF44'}
        
        fig = plt.figure(figsize=(14, 6))
        ax1 = fig.add_subplot(121, projection='3d')
        ax2 = fig.add_subplot(122)
        
        for i in range(n_atoms):
            col_offset = 1 + i * 4
            atom_id, atom_type = atom_info[i]
            x, y, z = data[:, col_offset+1], data[:, col_offset+2], data[:, col_offset+3]
            color = color_map.get(atom_type, '#888888')
            
            ax1.plot(x, y, z, color=color, linewidth=0.8, alpha=0.4, label=f'{atom_type} {atom_id}')
            ax2.plot(x, y, color=color, linewidth=0.8, alpha=0.4, label=f'{atom_type} {atom_id}')
        
        ax1.set_xlabel('X (Å)')
        ax1.set_ylabel('Y (Å)')
        ax1.set_zlabel('Z (Å)')
        ax1.set_title('3D Trajectories')
        if n_atoms <= 15:
            ax1.legend(fontsize=7)
        
        ax2.set_xlabel('X (Å)')
        ax2.set_ylabel('Y (Å)')
        ax2.set_title('XY Projection')
        ax2.set_aspect('equal')
        if n_atoms <= 15:
            ax2.legend(fontsize=7)
        
        plt.suptitle(f"Atom Trajectories - {self.sim_tag}")
        plt.tight_layout()
        self._save_plot("atom_trajectory")


class TrajectorySnapshot(AnalysisInput):
    """Trajectory snapshot generation."""
    FLAG = "trajectory_snapshot"
    CLI_HELP = "Generate trajectory snapshots"
    
    def __init__(self, data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps=None):
        super().__init__(data_dir, msd_dir, plots_dir, sim_tag, frame_interval_ps)
        self.CONFIG = ""
        self.filename = []
