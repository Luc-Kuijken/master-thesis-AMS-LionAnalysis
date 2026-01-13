#!/usr/bin/env amspython
# =====================================================================
# trajectory.py
# ---------------------------------------------------------------------
# Trajectory analysis runner using LionAnalysis with integrated plotting.
# Handles config generation, analysis execution, and visualization.
#
# Usage:
#   amspython trajectory.py TRAJ_FILE [options]
#
# Author: L. Kuijken
# Last updated: 2025-11-25
# =====================================================================

import argparse
import subprocess
import re
from datetime import datetime
from pathlib import Path

from analysis_inputs import AnalysisRegistry, AnalysisInput


# =====================================================================
# UTILITY FUNCTIONS
# =====================================================================
def log(message):
    """Print message with timestamp."""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {message}", flush=True)


# =====================================================================
# CLASS: Trajectory
# ---------------------------------------------------------------------
# Trajectory analysis runner for LionAnalysis.
# =====================================================================
class Trajectory:
    """
    Coordinates trajectory analysis workflow.
    
    - Generates LionAnalysis configuration files
    - Executes analysis
    - Generates plots using analysis classes
    """

    def __init__(
        self,
        traj_file,
        analysis_stride=1,
        n_threads=8,
        overwrite="skip",
    ):
        """
        Initialize Trajectory analyzer.
        
        Args:
            traj_file: Path to trajectory file
            analysis_stride: Stride for analysis (every Nth frame)
            n_threads: Number of threads for LionAnalysis
            overwrite: Overwrite mode ('skip', 'prompt', 'force')
        """
        self.traj_file = Path(traj_file)
        self.folder = self.traj_file.parent
        self.analysis_stride = analysis_stride
        self.n_threads = n_threads
        self.overwrite = overwrite

        # Parse simulation tag from filename
        self.sim_tag = self._parse_sim_tag()

        # Setup output directories
        self._setup_directories()

        # Load trajectory metadata
        self.metadata = self._load_metadata()

        # Length and timestep variables
        self.n_frames = int(self.metadata.get("n_frames", 10001))
        timestep_fs = float(self.metadata.get("output_sampling_fs", 100))
        self.timestep_ps = timestep_fs / 1000
        self.frame_interval_ps = self.timestep_ps * self.analysis_stride
        self.max_history = self.n_frames * self.timestep_ps

        # Regions for analysis (default to 'All Regions' if not specified)
        self.regions = self.metadata.get("regions", "All Regions")

        # Solvent indices (0-based, space-separated) - excludes MOF atoms if present
        self.solvent_indices = self._get_solvent_indices()

        # Initialize analysis registry
        self.analysis_registry = AnalysisRegistry(
            self.data_dir, 
            self.msd_dir, 
            self.plots_dir, 
            self.sim_tag, 
            frame_interval_ps=self.frame_interval_ps,
        )

    def __repr__(self):
        return (
            f"Trajectory(sim_tag={self.sim_tag}, "
            f"stride={self.analysis_stride}, "
            f"threads={self.n_threads}, "
            f"overwrite={self.overwrite})"
        )

    def _parse_sim_tag(self):
        """Extract simulation tag from trajectory filename."""
        stem = self.traj_file.stem
        if stem.startswith("MD_"):
            stem = stem[3:]
        if stem.endswith("_traj"):
            stem = stem[:-5]
        return stem

    def _setup_directories(self):
        """Create output directory structure."""
        base = Path.home() / "master_thesis_project" / "LionAnalysis"
        self.analysis_dir = base / self.folder.name
        self.data_dir = self.analysis_dir / "data"
        self.msd_dir = self.analysis_dir / "msd"
        self.plots_dir = self.analysis_dir / "plots"
        self.inout_dir = self.analysis_dir / "InOut"
        for directory in [
            self.analysis_dir,
            self.data_dir,
            self.msd_dir,
            self.plots_dir,
            self.inout_dir,
        ]:
            directory.mkdir(parents=True, exist_ok=True)

    def _load_metadata(self):
        """Load trajectory metadata from info file."""
        info_file = self.traj_file.with_name(f"{self.traj_file.stem}_xyz_info.txt")
        metadata = {}

        if info_file.exists():
            with open(info_file) as f:
                for line in f:
                    if ":" in line:
                        key, value = [x.strip() for x in line.split(":", 1)]
                        try:
                            metadata[key] = float(value)
                        except ValueError:
                            metadata[key] = value
        return metadata

    def _get_solvent_indices(self):
        """Get solvent atom indices (1-based, space-separated).
        """
        
        # Read total atom count from trajectory xyz file
        with open(self.traj_file) as f:
            total_atoms = int(f.readline().strip())
        
        # Parse H2O and H3O counts from sim_tag (e.g., "29H2O_4H3O_300K_100fs")
        h2o_match = re.search(r'(\d+)H2O', self.sim_tag)
        h3o_match = re.search(r'(\d+)H3O', self.sim_tag)
        
        n_h2o = int(h2o_match.group(1)) if h2o_match else 0
        n_h3o = int(h3o_match.group(1)) if h3o_match else 0
        
        # Calculate water atoms (H2O = 3 atoms, H3O = 4 atoms)
        n_water_atoms = (n_h2o * 3) + (n_h3o * 4)
        n_mof_atoms = total_atoms - n_water_atoms
        
        # Check if regions includes MOF
        regions = self.regions
        has_mof = "All Regions" in str(regions) or "MOF" in str(regions)
        
        if not has_mof or n_mof_atoms <= 0:
            # No MOF: all atoms are solvent (indices 1 to total)
            solvent_indices = " ".join(str(i) for i in range(1, total_atoms + 1))
            log(f"[Trajectory] No MOF: {total_atoms} solvent atoms (indices 1-{total_atoms})")
        else:
            # MOF present: solvent atoms start after MOF (indices n_mof+1 to total)
            solvent_indices = " ".join(str(i) for i in range(n_mof_atoms + 1, total_atoms + 1))
            log(f"[Trajectory] MOF detected: {n_mof_atoms} MOF atoms, {n_water_atoms} solvent atoms (indices {n_mof_atoms + 1}-{total_atoms})")
        
        return solvent_indices

    # -----------------------------------------------------------------
    # ANALYSIS STATUS CHECKING
    # -----------------------------------------------------------------
    def _get_expected_outputs(self, flag):
        """Get list of expected output files for an analysis type."""
        analysis = self.analysis_registry.get(flag)
        if not analysis:
            return []
        patterns = analysis.get_output_patterns()
        return [self.analysis_dir / p.format(tag=self.sim_tag) for p in patterns]

    def _is_analysis_complete(self, flag):
        """Check if all expected outputs exist for an analysis type."""
        expected = self._get_expected_outputs(flag)
        if not expected:
            return False

        missing = [p for p in expected if not p.exists()]
        if missing:
            log(f"[Trajectory] Missing outputs for '{flag}':")
            for p in missing:
                log(f"[Trajectory]    - {p.name}")
            return False
        return True

    # -----------------------------------------------------------------
    # CONFIG FILE GENERATION
    # -----------------------------------------------------------------
    def _generate_config(self, flags, atom_indices):
        """Generate LionAnalysis configuration file from analysis classes."""
        config_path = self.inout_dir / f"trajectory_{self.sim_tag}.config"

        # Handle atom indices: use provided value or default to "1" (valid placeholder for LionAnalysis)
        atom_indices_value = atom_indices if atom_indices else "1"

        # Base template variables
        template_vars = {
            "xyz_file": str(self.traj_file),
            "sim_tag": self.sim_tag,
            "output_dir": str(self.analysis_dir),
            "stride_length": self.analysis_stride,
            "timestep_ps": self.timestep_ps,
            "n_frames": self.n_frames,
            "max_history": self.max_history,
            "n_threads": self.n_threads,
            "solvent_indices": self.solvent_indices,
            "atom_indices": atom_indices_value,
            "stride": self.analysis_stride,
        }

        # Write combined config file
        config_parts = []
        
        # Add common config
        common_config = AnalysisInput.COMMON_CONFIG.format(**template_vars)
        config_parts.append(common_config)
        
        # Add group definitions (with atom_indices for tracking)
        group_def_config = AnalysisInput.GROUP_DEF_CONFIG.format(**template_vars)
        config_parts.append(group_def_config)
        
        # Add analysis-specific configs
        for flag_name, enabled in flags.items():
            if enabled:
                analysis = self.analysis_registry.get(flag_name)
                if analysis and analysis.CONFIG.strip():
                    # Create filename variables for this analysis (filename_0, filename_1, etc.)
                    for i, fname in enumerate(analysis.filename):
                        full_path = f"{analysis.output_dir}{fname}"
                        template_vars[f"filename_{i}"] = full_path
                    
                    # Format the CONFIG string with template variables
                    try:
                        formatted_analysis_config = analysis.CONFIG.format(**template_vars)
                        config_parts.append(formatted_analysis_config)
                    except KeyError as e:
                        log(f"[Trajectory] Warning: Missing variable {e} for {flag_name}")

        # Combine all parts
        full_config = "\n".join(config_parts)

        # Write to file
        with open(config_path, "w") as f_out:
            f_out.write(full_config)

        log(f"[Trajectory] Config generated for: {[f for f in flags if flags[f]]}")
        
        # Verify config was written
        if config_path.stat().st_size == 0:
            log(f"[Trajectory] ERROR: Generated config is empty!")
            return None
        
        return config_path

    # -----------------------------------------------------------------
    # LIONANALYSIS EXECUTION
    # -----------------------------------------------------------------
    def _run_lionanalysis(self, config_file, out_file):
        """Execute LionAnalysis with the generated config."""
        log(f"[Trajectory] Running LionAnalysis...")

        cmd = ["lionanalysis.exe", str(config_file)]
        with open(out_file, "w", encoding="utf-8") as fout:
            result = subprocess.run(
                cmd, stdout=fout, stderr=subprocess.STDOUT, text=True
            )

        if result.returncode != 0:
            log(f"[Trajectory] LionAnalysis failed (exit {result.returncode})")
            # Show last few lines of output for debugging
            try:
                with open(out_file, "r", encoding="utf-8") as f:
                    lines = f.read().splitlines()
                    for line in lines[-5:]:
                        log(f"  {line}")
            except Exception:
                pass
            raise subprocess.CalledProcessError(result.returncode, cmd)

    # -----------------------------------------------------------------
    # MAIN ANALYSIS RUNNER
    # -----------------------------------------------------------------
    def run(self, **flags):
        """Run trajectory analysis with specified flags."""
        # Extract atom_xyz from flags (not a boolean flag)
        atom_indices = flags.pop('atom_xyz', None)
        
        log(f"[Trajectory]")
        log(f"[Trajectory] Analysis: {self.sim_tag}")
        log(f"[Trajectory] Output: {self.analysis_dir}")
        log(f"[Trajectory] Flags: {[f for f in flags if flags[f]]}")
        if atom_indices:
            log(f"[Trajectory] Tracking atoms: {atom_indices}")

        # Determine which analyses to run
        active_flags = self._get_active_analyses(flags)

        # Generate plots if no analysis needed
        if not active_flags:
            log(f"[Trajectory] No new analysis needed")
            self._generate_plots(flags)
            return

        # Generate config and run analysis
        config_path = self._generate_config(active_flags, atom_indices)
        if not config_path:
            log("[Trajectory] Config generation failed")
            return

        out_file = self.inout_dir / f"trajectory_{self.sim_tag}.out"
        self._run_lionanalysis(config_path, out_file)
        if config_path.exists():
            config_path.unlink()  # Clean up config file

        log("[Trajectory] Analysis complete")

        # Verify outputs and generate plots
        self._verify_outputs(active_flags)
        self._generate_plots(flags)

    def _get_active_analyses(self, flags):
        """Determine which analyses need to run based on completion status."""
        active_flags = {}

        for key, enabled in flags.items():
            if not enabled:
                continue

            if self._is_analysis_complete(key) and self.overwrite != "force":
                if self.overwrite == "skip":
                    log(f"[Trajectory] Skipping '{key}': already complete")
                    continue
                elif self.overwrite == "prompt":
                    response = input(
                        f"'{key}' outputs exist. Re-run? [y/N] "
                    ).strip().lower()
                    if response not in {"y", "yes"}:
                        continue

            active_flags[key] = enabled

        return active_flags

    def _verify_outputs(self, flags):
        """Verify that expected outputs were created and run post-processing."""
        for flag in flags:
            # Run post-processing if defined (e.g., rename files)
            analysis = self.analysis_registry.get(flag)
            if analysis and hasattr(analysis, 'post_process'):
                analysis.post_process()
            
            # Check if outputs are complete
            if not self._is_analysis_complete(flag):
                log(f"[Trajectory] Warning: Some '{flag}' outputs missing")

    def _generate_plots(self, flags):
        """Generate plots from analysis outputs using analysis classes."""
        try:
            # Call plot method for each enabled analysis
            for flag_name, enabled in flags.items():
                if enabled:
                    analysis = self.analysis_registry.get(flag_name)
                    if analysis:
                        try:
                            log(f"[Plotter] Generating plots for {flag_name}...")
                            analysis.plot()
                        except Exception as e:
                            log(f"[Plotter] Failed to plot {flag_name}: {e}")
            
            log(f"[Plotter] Finished plotting for {self.sim_tag}")
        except Exception as e:
            log(f"[Trajectory] Plotting failed: {e}")


# =====================================================================
# CLI INTERFACE
# =====================================================================
def build_parser():
    """Build command-line argument parser from registered analysis classes."""
    parser = argparse.ArgumentParser(
        description="Run LionAnalysis trajectory analysis with integrated plotting."
    )
    parser.add_argument("traj_file", type=str, help="Path to trajectory file")
    
    # Auto-add arguments from registered analysis classes
    for flag, cls in AnalysisInput.get_all_analysis_classes().items():
        if cls.CLI_ACTION == "store_true":
            parser.add_argument(
                f"--{flag}", 
                action="store_true", 
                help=cls.CLI_HELP
            )
        elif cls.CLI_ACTION == "store":
            parser.add_argument(
                f"--{flag}", 
                type=str, 
                metavar=cls.CLI_VAR, 
                help=cls.CLI_HELP
            )
    
    # Common arguments
    parser.add_argument(
        "--overwrite",
        choices=["prompt", "skip", "force"],
        default="skip",
        help="Overwrite mode for existing outputs",
    )
    parser.add_argument(
        "--stride", 
        type=int, 
        default=1, 
        help="Analysis stride"
    )
    parser.add_argument(
        "--n_threads", 
        type=int, 
        default=8, 
        help="Number of threads"
    )
    return parser


def main(argv=None):
    """Main entry point for command-line usage."""
    args = build_parser().parse_args(argv)

    traj = Trajectory(
        args.traj_file,
        args.stride,
        args.n_threads,
        args.overwrite,
    )

    # Build flags dictionary from registered analysis classes
    flags = {}
    for flag, cls in AnalysisInput.get_all_analysis_classes().items():
        arg_value = getattr(args, flag, None)
        if cls.CLI_ACTION == "store_true":
            flags[flag] = arg_value if arg_value else False
        elif cls.CLI_ACTION == "store":
            if arg_value:
                flags[flag] = ' '.join(arg_value.replace(',', ' ').split())
            else:
                flags[flag] = None

    traj.run(**flags)


if __name__ == "__main__":
    main()