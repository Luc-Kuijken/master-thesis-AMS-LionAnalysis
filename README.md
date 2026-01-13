# Master Thesis Project: Trajectory Analysis for Molecular Dynamics Simulations

A comprehensive toolkit for analyzing molecular dynamics trajectories of water and hydronium systems in Metal-Organic Frameworks (MOFs), specifically UiO-66.

## Project Overview

This project provides tools for:
- **Trajectory analysis**: RDF, MSD, hydrogen bonding, proton transfer events
- **Parallel processing**: Efficient batch analysis of multiple trajectories
- **Visualization**: Automated plotting and data export
- **Simulation management**: AMS/SCM workflow integration

## Project Structure

```
master_thesis_project/
├── Tools/                      # Core analysis modules
│   ├── trajectory.py           # Main trajectory analysis script
│   ├── analysis_inputs.py      # Configuration and input handling
│   ├── common.py               # Shared utilities
│   ├── simulation.py           # Simulation setup and management
│   └── configs/                # Analysis configuration files
│
├── AMS_simulations/            # Molecular dynamics simulation data
│   ├── MD_33H2O_0H3O/         # Small system (33 water molecules)
│   ├── MD_65H2O_0H3O/         # Medium system (65 water molecules)
│   ├── MD_98H2O_0H3O/         # Large system (98 water molecules)
│   └── MD_130H2O_0H3O/        # XL system (130 water molecules)
│
├── Trajectories/               # Converted trajectory files
│   ├── xyz/                    # XYZ format trajectories
│   └── pdb/                    # PDB format trajectories
│
├── LionAnalysis/               # Analysis results (Lion = Li+ ion studies)
│   ├── MD_33H2O/              # Results for small system
│   ├── MD_65H2O/              # Results for medium system
│   ├── MD_98H2O/              # Results for large system
│   └── MD_130H2O/             # Results for XL system
│
└── SitesAnalysis/              # Site-specific analysis results
    ├── AvOPs/                  # Average occupancy positions
    └── data_files/             # Processed data
```

## Key Features

### 1. Trajectory Analysis (`trajectory.py`)
Analyzes molecular dynamics trajectories with support for:
- **Radial Distribution Functions (RDF)**: Pair correlations
- **Mean Square Displacement (MSD)**: Water and hydronium diffusion
- **Hydrogen Bond Analysis**: Network dynamics
- **Proton Transfer Events**: PTFEL (Proton Transfer Free Energy Landscape)
- **Residence Times**: Ion-site interactions
- **Atom Tracking**: Individual particle trajectories

### 2. Parallel Processing
The `zanalyse_trajectories.sh` script enables:
- Batch processing of multiple trajectories
- Configurable parallelism (jobs and threads)
- Progress tracking and logging
- Error recovery

### 3. System Configurations
Supports multiple system sizes and conditions:
- **Water content**: 33, 65, 98, 130 H₂O molecules
- **Hydronium ions**: 0-4 H₃O⁺ ions
- **Temperatures**: 300K, 350K, 400K
- **Time steps**: 20fs, 100fs

## Installation

### Prerequisites
- Python 3.8+
- AMS/SCM software suite (for `amspython`)
- NumPy, SciPy, Matplotlib

### Setup
```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/master-thesis-trajectory-analysis.git
cd master-thesis-trajectory-analysis

# The project uses amspython (part of AMS/SCM suite)
# Ensure AMS is installed and sourced:
source /path/to/ams/amsbashrc.sh
```

## Usage

### Basic Trajectory Analysis
```bash
# Analyze a single trajectory
amspython Tools/trajectory.py Trajectories/xyz/MD_130H2O/MD_130H2O_0H3O_300K_100fs_traj.xyz --rdf --msd_water

# Run with all analysis options
amspython Tools/trajectory.py trajectory.xyz --rdf --msd_water --msd_hydronium --h_bonds --ptfel --res_times
```

### Batch Processing
```bash
# Analyze all trajectories in a directory (parallel processing)
zanalyse_trajectories.sh Trajectories/xyz/MD_130H2O/ --rdf --msd_water --threads 8 --jobs 40

# Force overwrite existing results
zanalyse_trajectories.sh Trajectories/xyz/ --rdf --overwrite force --threads 4 --jobs 60
```

### Analysis Options
```
--rdf                Enable radial distribution function analysis
--msd_water          Enable water mean square displacement
--msd_hydronium      Enable hydronium mean square displacement
--h_bonds            Enable hydrogen bond analysis
--ptfel              Enable proton transfer analysis
--res_times          Enable residence time analysis
--atom_xyz INDEX     Track specific atom by index
--stride N           Analysis stride (frames per step)
--threads N          Threads per job
--overwrite MODE     Overwrite mode: prompt, skip, or force
```

## Output

Analysis results are organized by system and analysis type:
```
LionAnalysis/MD_130H2O/
├── rdf_data/                   # Radial distribution functions
├── msd_data/                   # Mean square displacement
├── hbond_data/                 # Hydrogen bond statistics
├── ptfel_data/                 # Proton transfer events
└── plots/                      # Generated visualizations
```

## Configuration

Analysis parameters are managed through:
- `Tools/configs/`: YAML/JSON configuration files
- `Tools/analysis_inputs.py`: Input validation and parsing
- Command-line arguments: Override defaults per run

## Logging

Comprehensive logging system:
- **Master log**: `LionAnalysis/Analysis_logs/master_summary.log`
- **Job logs**: Individual trajectory analysis logs
- **Error tracking**: Failed analyses recorded for retry

## Contributing

This is a research project. For questions or collaboration:
- Contact: [Your Name/Email]
- Institution: [Your Institution]

## License

[Specify your license - e.g., MIT, GPL, or Academic/Research Use Only]

## Acknowledgments

- AMS/SCM software suite for molecular dynamics simulations
- [Your advisor/research group]
- [Funding sources]

## Citation

If you use this code in your research, please cite:
```
[Your thesis citation to be added]
```
