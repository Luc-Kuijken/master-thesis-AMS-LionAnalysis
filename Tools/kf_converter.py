#!/usr/bin/env amspython
# =====================================================================
# kf_converter.py
# ---------------------------------------------------------------------
# Convert AMS RKF trajectories into XYZ or PDB trajectories.
# Supports region filtering, sampling adjustment, and max-time settings.
#
# Usage:
#   amspython kf_converter.py FOLDER [options]
#
# Author: L. Kuijken
# Last updated: 2025-11-26
# =====================================================================

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
from typing import Iterable, List, Optional

import sys

from scm.plams import KFFile, PeriodicTable, Trajectory as PlamsTrajectory, Units


# =====================================================================
# UTILITY FUNCTIONS
# =====================================================================
def log(message: str) -> None:
    """Print message with timestamp"""
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {message}", flush=True)


def parse_time_str(time_str: Optional[str]) -> Optional[int]:
    """Convert time strings like '100fs', '2ps', '0.5ns' to integer fs"""
    if not time_str:
        return None
    val = time_str.strip().lower()
    if val.endswith("fs"):
        return int(float(val[:-2]))
    elif val.endswith("ps"):
        return int(float(val[:-2]) * 1e3)
    elif val.endswith("ns"):
        return int(float(val[:-2]) * 1e6)
    raise ValueError("Time inputs must include units (fs, ps, ns).")


def parse_region_list(region: Optional[str]) -> Optional[List[str]]:
    """Parse comma-separated region list string (moved from old _parse_region)."""
    if not region:
        return None
    return [r.strip() for r in region.split(",") if r.strip()]


# =====================================================================
# CLASS: KFReader
# ---------------------------------------------------------------------
# Handles reading from the AMS RKF trajectory file.
# =====================================================================
class KFReader:
    def __init__(self, kf_path: Path):
        self.kf_path = kf_path
        self.reader = KFFile(str(kf_path))
        self.bohr_to_ang = Units.conversion_ratio("bohr", "angstrom")

    def get_sampling_interval(self) -> int:
        """Return the sampling interval (fs) from the RKF file."""
        step1 = self.reader.read("History", "Step(1)")
        step2 = self.reader.read("History", "Step(2)")
        return int(step2 - step1)

    def compute_md_timestep(self) -> Optional[float]:
        """Compute the MD timestep from MDResults, if available."""
        try:
            start_time = self.reader.read("MDResults", "StartTime[fs]")
            end_time = self.reader.read("MDResults", "EndTime[fs]")
            start_step = self.reader.read("MDResults", "StartStep")
            end_step = self.reader.read("MDResults", "EndStep")
        except Exception:
            return None
        step_span = end_step - start_step
        if step_span <= 0:
            return None
        return (end_time - start_time) / float(step_span)

    def get_atom_data(self):
        """Return atomic symbols and lattice vectors."""
        atom_numbers = self.reader.read("Molecule", "AtomicNumbers")
        symbols = [PeriodicTable.get_symbol(n) for n in atom_numbers]
        lattice = [v * self.bohr_to_ang for v in self.reader.read("Molecule", "LatticeVectors")]
        return symbols, lattice

    def get_coords(self, frame: int):
        """Return coordinates for the given frame index."""
        return self.reader.read("History", f"Coords({frame})")
    

# =====================================================================
# CLASS: TrajectoryWriter
# ---------------------------------------------------------------------
# Handles writing XYZ/PDB trajectory frames and metadata files.
# =====================================================================
class TrajectoryWriter:
    def __init__(self, out_base: Path, fmt: str, sim_folder_name: str):
        self.out_base = out_base
        self.fmt = fmt
        self.sim_folder_name = sim_folder_name
        self._last_progress = None
    
    def progress(self, frame: int, total: int) -> None:
        """Display conversion progress every 10%."""
        progress = int((frame / total) * 100)
        if progress % 10 == 0 and progress > 0:
            if self._last_progress != progress:
                log(f"[KFConverter] {self.sim_folder_name}: {progress}% complete ({frame}/{total})")
                self._last_progress = progress

    # -----------------------------------------------------------------
    # TRAJECTORY WRITERS (Additional methods for writing frames can be added here)
    # -----------------------------------------------------------------
    def write_xyz_frames(self, reader: KFReader, indices, symbols, lattice, unit_conversion, stride, n_frames):
        """Write trajectory frames in XYZ format."""
        frames: List[str] = []
        for frame in range(1, n_frames + 1, stride):
            coords = reader.get_coords(frame)
            a, b, c = lattice[0], lattice[4], lattice[8]
            frame_data = [f"{len(indices)}", f"XYZ {a} {b} {c}"]
            for i in indices:
                x, y, z = (coords[3 * i + j] * unit_conversion for j in range(3))
                frame_data.append(f"{symbols[i]} {x:8.5f} {y:8.5f} {z:8.5f}")
            frames.append("\n".join(frame_data))
            self.progress(frame, n_frames)
        return frames

    def write_pdb_frames(self, reader: KFReader, indices, symbols, lattice, unit_conversion, stride, n_frames):
        """Write trajectory frames in PDB format."""
        frames: List[str] = []
        for frame in range(1, n_frames + 1, stride):
            coords = reader.get_coords(frame)
            a, b, c = lattice[0], lattice[4], lattice[8]
            frame_data = [f"MODEL   {frame}", f"CRYST1 {a:8} {b:8} {c:8} 90.00 90.00 90.00"]
            for idx, atom_index in enumerate(indices):
                x, y, z = (coords[3 * atom_index + j] * unit_conversion for j in range(3))
                frame_data.append(
                    f"ATOM {idx + 1:6d} {symbols[atom_index]:<4} MOL     {x:8.3f}{y:8.3f}{z:8.3f}"
                )
            frame_data.append("ENDMDL")
            frames.append("\n".join(frame_data))
            self.progress(frame, n_frames)
        return frames

    def _reorder_h2o_for_pdb(self, indices, symbols, regions) -> List[int]:
        """Reorder H2O triplets from H-O-H to O-H-H pattern for PDB output."""
        if not indices:
            return []
        reordered: List[int] = []
        i = 0
        while i < len(indices):
            if (
                i + 2 < len(indices)
                and regions[indices[i]] == "H2O"
                and regions[indices[i + 1]] == "H2O"
                and regions[indices[i + 2]] == "H2O"
                and symbols[indices[i]] == "H"
                and symbols[indices[i + 1]] == "O"
                and symbols[indices[i + 2]] == "H"
            ):
                reordered.extend([indices[i + 1], indices[i], indices[i + 2]])
                i += 3
            else:
                reordered.append(indices[i])
                i += 1
        return reordered

    def write_info(
        self,
        info_path: Path,
        md_timestep_fs: float,
        saved_sampling_fs: float,
        eff_sampling_fs: float,
        n_frames: int,
        total_time_fs: float,
        regions: List[str],
    ) -> None:
        """Write metadata file for the trajectory (identical output as original)."""
        with open(info_path, "w") as handle:
            handle.write(f"md_timestep_fs: {int(round(md_timestep_fs))}\n")
            handle.write(f"saved_sampling_fs: {saved_sampling_fs}\n")
            handle.write(f"output_sampling_fs: {eff_sampling_fs}\n")
            handle.write(f"n_frames: {n_frames}\n")
            handle.write(f"simulation_length_fs: {int(round(total_time_fs))}\n")
            handle.write(f"regions: {regions}\n")
    

# =====================================================================
# CLASS: KFConverter
# ---------------------------------------------------------------------
# Handles reading, filtering, conversion, and writing.
# =====================================================================
class KFConverter:
    def __init__(self, folder_path: Path, out_base=None, fmt="matti", overwrite_mode="prompt"):
        self.sim_folder = Path(folder_path).expanduser().resolve()
        self._get_folder_metadata()

        self.kf_file = self.sim_folder / "ams.results" / "ams.rkf"
        if not self.kf_file.exists():
            raise FileNotFoundError(f"RKF not found: {self.kf_file}")

        self.fmt = fmt.lower()
        default_base = f"~/master_thesis_project/Trajectories/{'xyz' if self.fmt == 'matti' else 'pdb'}"
        self.out_base = Path(out_base or default_base).expanduser().resolve()
        self.out_base.mkdir(parents=True, exist_ok=True)

        if overwrite_mode not in {"prompt", "skip", "force"}:
            raise ValueError("overwrite_mode must be 'prompt', 'skip', or 'force'")
        self.overwrite_mode = overwrite_mode

    def __repr__(self):
        return (
            f"KFConverter(Total molecules={self.n_total}, H2O={self.n_H2O}, "
            f"H3O={self.n_H3O}, Temp={self.temperature}, "
            f"OrigFreq={self.original_freq}fs, Format={self.fmt})"
        )

    def _get_folder_metadata(self) -> None:
        """Extract H2O/H3O counts, temperature, and sampling frequency from folder name."""
        parts = self.sim_folder.name.split("_")
        self.n_H2O = int(parts[1].replace("H2O", ""))
        self.n_H3O = int(parts[2].replace("H3O", ""))
        self.n_total = self.n_H2O + self.n_H3O
        self.temperature = int(parts[3].replace("K", ""))
        self.original_freq = parse_time_str(parts[4])

    def _region_indices(self, traj: PlamsTrajectory, region: Optional[List[str]]):
        """Return indices of atoms matching region filter."""
        if not region:
            return range(len(traj[0]))
        return [
            i for i, atom in enumerate(traj[0])
            if hasattr(atom.properties, "region") and atom.properties.region in region
        ]

    # -----------------------------------------------------------------
    # MAIN CONVERSION FUNCTION
    # -----------------------------------------------------------------
    def convert(
        self,
        region: Optional[str] = None,
        output_sampling: Optional[str] = None,
        max_time: Optional[str] = None,
        fmt: Optional[str] = None,
        overwrite_mode: Optional[str] = None,
    ) -> Path:
        """Convert RKF trajectory to specified format with filtering options."""
        region_list = parse_region_list(region)
        output_sampling_fs = parse_time_str(output_sampling) if output_sampling else None
        max_time_fs = parse_time_str(max_time) if max_time else None
        fmt = (fmt or self.fmt).lower()
        mode = overwrite_mode or self.overwrite_mode

        # Automatically switch to skip mode in non-interactive environments (e.g. SLURM)
        if mode == "prompt" and not sys.stdin.isatty():
            mode = "skip"
            log(f"[KFConverter] Non-interactive environment detected -> switching overwrite mode to 'skip'")

        # Initialize readers and writers
        reader = KFReader(self.kf_file)
        traj = PlamsTrajectory(str(self.kf_file))
        writer = TrajectoryWriter(self.out_base, fmt, self.sim_folder.name)

        sampling_fs = reader.get_sampling_interval()
        md_timestep_fs = reader.compute_md_timestep() or sampling_fs

        # Determine stride and new output frames
        stride = self._determine_stride(sampling_fs, output_sampling_fs)
        total_frames = len(traj)
        if max_time_fs:
            max_frame = int(max_time_fs // sampling_fs) + 1
            total_frames = min(total_frames, max_frame)
        output_frames = ((total_frames - 1) // stride) + 1 if total_frames else 0
        last_frame_index = 1 + (output_frames - 1) * stride if output_frames else 0

        # Check if files already exist before processing
        out_folder, traj_file, info_file = self._prepare_output(fmt, output_sampling_fs)
        if self._should_skip(mode, traj_file, info_file, out_folder):
            return traj_file  # Return path to existing trajectory file

        symbols, lattice = reader.get_atom_data()
        region_indices = list(self._region_indices(traj, region_list)) or list(range(len(traj[0])))
        region_labels = [getattr(atom.properties, "region", None) for atom in traj[0]]

        # Write frames
        if fmt == "matti":
            frames = writer.write_xyz_frames(reader, region_indices, symbols, lattice, reader.bohr_to_ang, stride, last_frame_index)
        elif fmt == "pdb":
            reordered = writer._reorder_h2o_for_pdb(region_indices, symbols, region_labels)
            frames = writer.write_pdb_frames(reader, reordered or region_indices, symbols, lattice, reader.bohr_to_ang, stride, last_frame_index)
        else:
            raise ValueError(f"Unknown format: {fmt}")

        with open(traj_file, "w") as f:
            f.write("\n".join(frames))

        total_time_fs = (last_frame_index - 1) * sampling_fs if output_frames else 0.0
        writer.write_info(info_file, md_timestep_fs, sampling_fs, stride * sampling_fs, output_frames, total_time_fs, region_list or ["All Regions"])

        log(f"[KFConverter] Wrote {traj_file.name} ({output_frames} frames, stride={stride})")
        return traj_file  # Return path to trajectory file

    # -----------------------------------------------------------------
    # HELPER METHODS
    # -----------------------------------------------------------------
    def _determine_stride(self, sampling_fs: int, output_sampling_fs: Optional[int]) -> int:
        """Determine stride factor from requested output sampling."""
        if not output_sampling_fs:
            return 1
        ratio = output_sampling_fs / sampling_fs
        stride = max(int(round(ratio)), 1)
        if abs(ratio - stride) > 1e-6:
            raise ValueError(f"Requested output sampling {output_sampling_fs} fs not multiple of {sampling_fs} fs")
        if ratio < 1 - 1e-6:
            raise ValueError(f"Requested output sampling {output_sampling_fs} fs finer than {sampling_fs} fs")
        return stride

    def _prepare_output(self, fmt: str, output_sampling_fs: Optional[int]) -> tuple[Path, Path, Path]:
        """Prepare output file paths."""
        out_folder = self.out_base / f"MD_{self.n_total}H2O"
        out_folder.mkdir(parents=True, exist_ok=True)
        ext = "xyz" if fmt == "matti" else "pdb"
        traj_file = out_folder / f"MD_{self.n_H2O}H2O_{self.n_H3O}H3O_{self.temperature}K_{output_sampling_fs or self.original_freq}fs_traj.{ext}"
        info_file = traj_file.with_name(f"{traj_file.stem}_{ext}_info.txt")
        return out_folder, traj_file, info_file

    def _should_skip(self, mode: str, traj_file: Path, info_file: Path, out_folder: Path) -> bool:
        """Handle overwrite decisions."""
        if not (traj_file.exists() or info_file.exists()):
            return False
        if mode == "prompt":
            ans = input(f"{traj_file.name} exists in {out_folder}. Overwrite? [y/N] ").strip().lower()
            mode = "force" if ans in {"y", "yes"} else "skip"
        if mode == "skip":
            log(f"[KFConverter] Skipped conversion for {self.sim_folder.name}")
            return True
        return False

# =====================================================================
# CLI INTERFACE
# =====================================================================
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert AMS RKF trajectory to XYZ or PDB format.")
    parser.add_argument("folder_path", type=str, help="Path to AMS simulation folder.")
    parser.add_argument("-r", "--region", type=str, help="Region(s) to include (comma-separated).")
    parser.add_argument("-s", "--output_sampling", type=str, help="Output sampling interval (fs, ps, ns).")
    parser.add_argument("-t", "--max_time", type=str, help="Maximum simulation time (fs, ps, ns).")
    parser.add_argument("-f", "--format", dest="fmt", type=str, default="matti", choices=["matti", "pdb"])
    parser.add_argument("-o", "--out", type=str, help="Output base directory.")
    parser.add_argument("--overwrite", dest="overwrite_mode", type=str, default="prompt", choices=["prompt", "skip", "force"])
    return parser


def main(argv: Optional[Iterable[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    converter = KFConverter(Path(args.folder_path), args.out, args.fmt, args.overwrite_mode)
    traj_file = converter.convert(region=args.region, output_sampling=args.output_sampling, max_time=args.max_time)
    log(f"[KFConverter] Trajectory saved at {traj_file}")


if __name__ == "__main__":
    main()
