#!/usr/bin/env amspython
# =====================================================================
# simulation.py
# ---------------------------------------------------------------------
# Creates AMS MD simulation folders with tagged molecule structures
# and run.run files. Supports MOPS, HPC, and Snellius configurations.
#
# Usage:
#   amspython simulation.py [options]
#
# Author: L. Kuijken
# Last updated: 2025-11-13
# =====================================================================

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path
from typing import Optional

import numpy as np
from scm.plams import AMSJob, Atom, Molecule, Settings, packmol_in_void


# =====================================================================
# UTILITY FUNCTIONS
# =====================================================================
def log(message: str) -> None:
    """Print message with timestamp."""
    from datetime import datetime
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {message}", flush=True)


# =====================================================================
# CLASS: Structure
# ---------------------------------------------------------------------
# Handles MOF structure creation, solvent packing, and functionalization.
# =====================================================================
class Structure:
    """Create and manipulate MOF structures with solvents."""

    def __init__(self, mof_file: str):
        self.mof_path = Path(mof_file)
        self.name = self.mof_path.stem
        self.output_file = None
        self.tagged_file = None
        self.idx = []
        self.region = []
        self.n_atoms = 0

        if not self.mof_path.is_file():
            raise FileNotFoundError(f"MOF file not found: {mof_file}")

        self.struct = Molecule(str(self.mof_path))
        self.idx = [len(self.struct)]
        self.region = ["MOF"]
        self.solvent_counts = []
        self.solvent_names = []
        self.read_unitcell_from_xyz(mof_file)

        log(f"Loaded MOF structure from {self.mof_path.stem} with {len(self.struct)} atoms.")

    def __repr__(self):
        return f"Structure(name={self.name}, n_atoms={len(self.struct)}, solvents={self.solvent_names})"

    def read_unitcell_from_xyz(self, filename):
        """Read unit cell vectors from AMS-style XYZ file."""
        vecs = []
        with open(filename, "r") as f:
            for line in f:
                if line.startswith("VEC"):
                    parts = line.split()
                    vec = [float(x) for x in parts[1:4]]
                    vecs.append(vec)
        if len(vecs) == 3:
            self.unitcell = np.array(vecs)
        else:
            raise ValueError("Could not read 3 unit cell vectors from xyz")

    def displacement_pbc(self, a, b, cell):
        """Return minimum-image displacement vector (a→b) under PBC."""
        diff = np.array(b.coords) - np.array(a.coords)
        frac = np.linalg.solve(cell.T, diff)
        frac -= np.round(frac)
        return cell.T.dot(frac)

    def close_atoms(self, atom, mol, cutoff=1.6):
        """Return all atoms within cutoff Å of atom (except itself), using PBC."""
        neighbors = []
        for n in mol:
            if n is atom:
                continue
            dist = np.linalg.norm(self.displacement_pbc(atom, n, self.unitcell))
            if dist < cutoff:
                neighbors.append(n)
        return neighbors

    def add_amino_group(self, mol, carbon):
        """Add NH2 group to aromatic carbon."""
        # Remove H bound to carbon
        H_neighbors = [h for h in self.close_atoms(carbon, mol, cutoff=1.2) if h.symbol == "H"]
        if len(H_neighbors) != 1:
            log(f"Warning: expected 1 H on C{mol.index(carbon)}, found {len(H_neighbors)}")
            return
        h_atom = H_neighbors[0]
        Hpos = np.array(h_atom.coords)
        mol.delete_atom(h_atom)

        # Get ring neighbors
        ring_neighbors = [n for n in self.close_atoms(carbon, mol, cutoff=1.6) if n.symbol == "C"]
        if len(ring_neighbors) < 2:
            log(f"Warning: expected >=2 ring neighbors on C{mol.index(carbon)}, found {len(ring_neighbors)}")
            return

        v1 = self.displacement_pbc(carbon, ring_neighbors[0], self.unitcell)
        v2 = self.displacement_pbc(carbon, ring_neighbors[1], self.unitcell)

        # Calculate ring plane normal
        normal = np.cross(v1, v2)
        normal /= np.linalg.norm(normal)

        # Place N and H atoms
        Npos = Hpos
        N = Atom(symbol="N", coords=Npos)
        mol.add_atom(N)
        mol.add_bond(carbon, N, order=1)

        H1 = Atom(symbol="H", coords=Npos + normal * 1.0)
        H2 = Atom(symbol="H", coords=Npos - normal * 1.0)
        mol.add_atom(H1)
        mol.add_atom(H2)
        mol.add_bond(N, H1, order=1)
        mol.add_bond(N, H2, order=1)

    def find_linkers(self, mol):
        """Identify all unique terephthalate linkers in UiO-66."""
        linkers = []
        seen = set()

        # Find carboxylate carbons
        carboxylates = []
        assigned_carboxylates = set()
        for a in mol:
            if a.symbol == "C":
                neighbors = self.close_atoms(a, mol, cutoff=1.6)
                O_close = [n for n in neighbors if n.symbol == "O"]
                if len(O_close) >= 2:
                    carboxylates.append(a)

        # Group into linkers
        for c1 in carboxylates:
            if c1 in assigned_carboxylates:
                continue

            for neigh in self.close_atoms(c1, mol, cutoff=1.6):
                if neigh.symbol == "C":
                    ring_atoms = set()
                    to_visit = [neigh]

                    while to_visit:
                        current = to_visit.pop()
                        if current in ring_atoms:
                            continue
                        ring_atoms.add(current)

                        for nn in self.close_atoms(current, mol, cutoff=1.6):
                            if nn.symbol == "C" and nn not in ring_atoms:
                                to_visit.append(nn)

                    carboxylates_in_ring = [c for c in ring_atoms if c in carboxylates]

                    if len(carboxylates_in_ring) == 2 and len(ring_atoms) >= 6:
                        linker_atoms = ring_atoms.union(carboxylates_in_ring)
                        sig = tuple(sorted([mol.index(a) for a in linker_atoms]))
                        if sig not in seen:
                            seen.add(sig)
                            linkers.append(linker_atoms)
                            assigned_carboxylates.update(carboxylates_in_ring)

        return linkers

    def find_ortho_carbons_per_linker(self):
        """For each linker, return one carbon at the ortho position."""
        linkers = self.find_linkers(self.struct)
        ortho_sites = []

        for linker in linkers:
            candidates = []
            for atom in linker:
                if atom.symbol == "C":
                    O_neighbors = [n for n in self.close_atoms(atom, self.struct, cutoff=1.6) if n.symbol == "O"]
                    if len(O_neighbors) >= 2:  # carboxylate carbon
                        for neigh in self.close_atoms(atom, self.struct, cutoff=1.6):
                            if neigh.symbol == "C" and neigh in linker:
                                for nn in self.close_atoms(neigh, self.struct, cutoff=1.6):
                                    if nn.symbol == "C" and nn in linker and nn is not atom:
                                        candidates.append(nn)
            if candidates:
                chosen = min(candidates, key=lambda x: self.struct.index(x))
                ortho_sites.append(chosen)

        return ortho_sites

    def functionalize_linkers(self, group="NH2"):
        """Functionalize linkers with specified group."""
        ortho_sites = self.find_ortho_carbons_per_linker()
        log(f"Found {len(ortho_sites)} ortho carbons to functionalize.")

        if group == "NH2":
            for carbon in ortho_sites:
                self.add_amino_group(self.struct, carbon)

        log(f"Functionalized {len(ortho_sites)} sites with {group}.")

    def pack_solvent(self, solvent_file, n_mol, tolerance: float = 2.0):
        """Pack solvent molecules into MOF structure."""
        solvent_path = Path(solvent_file)
        if not solvent_path.is_file():
            raise FileNotFoundError(f"Solvent file not found: {solvent_file}")

        solvent = Molecule(str(solvent_path))
        len_solvent = len(solvent)

        self.struct = packmol_in_void(
            host=self.struct,
            molecules=solvent,
            n_molecules=n_mol,
            tolerance=tolerance,
        )

        log(f"Packed {n_mol} molecules of {solvent_path.stem} into the structure.")

        self.idx.append(self.idx[-1] + n_mol * len_solvent)
        self.region.append(solvent_path.stem)
        self.solvent_counts.append(n_mol)
        self.solvent_names.append(solvent_path.stem)

    def name_output_file(self):
        """Generate output filename based on composition."""
        name_parts = [f"{n}_{s}" for n, s in zip(self.solvent_counts, self.solvent_names)]
        suffix = "_".join(name_parts)
        self.output_file = f"{self.name}_{suffix}.xyz"
        self.struct.write(self.output_file, outputformat="xyz")

    def write_regions_and_charge(self):
        """Write tagged XYZ file with region information."""
        if not self.output_file:
            raise ValueError("Solvent must be packed before tagging regions.")

        self.tagged_file = self.output_file.replace(".xyz", "_tagged.xyz")

        with open(self.output_file, "r") as f:
            lines = f.readlines()

        header = lines[:2]
        atom_lines = lines[2:-3]
        footer = lines[-3:]

        new_lines = []
        i = 0
        for j in range(len(self.idx)):
            while i < self.idx[j]:
                line = atom_lines[i].rstrip("\n")
                new_lines.append(f"{line} region={self.region[j]}\n")
                i += 1

        with open(self.tagged_file, "w") as f:
            f.writelines(header)
            f.writelines(new_lines)
            f.writelines(footer)

    def get_output_files(self):
        """Return paths to output files."""
        return self.output_file, self.tagged_file


# =====================================================================
# CLASS: Simulation
# ---------------------------------------------------------------------
# Handles simulation folder creation and run.run file generation.
# =====================================================================
class Simulation:
    """Create AMS MD simulation folders."""

    def __init__(
        self,
        base_path,
        n_H2O,
        n_H3O,
        temperature,
        mof_file,
        h2o_file,
        h3o_file,
        n_steps=10000,
        sampling_freq=100,
        sim_threads=32,
        NSCM=1,
        machine="MOPS",
    ):
        self.n_H2O = n_H2O
        self.n_H3O = n_H3O
        self.temperature = temperature
        self.base_path = Path(base_path)
        self.mof_file = mof_file
        self.h2o_file = h2o_file
        self.h3o_file = h3o_file

        self.n_steps = n_steps
        self.sampling_freq = sampling_freq
        self.sim_threads = sim_threads
        self.NSCM = NSCM
        self.machine = machine.upper()

        self.folder_path = self.base_path / f"MD_{self.n_H2O}H2O_{self.n_H3O}H3O" / f"MD_{self.n_H2O}H2O_{self.n_H3O}H3O_{self.temperature}K_{self.sampling_freq}fs"
        self.name = f"MD_{self.n_H2O}H2O_{self.n_H3O}H3O_{self.temperature}K_{self.sampling_freq}fs"

    def __repr__(self):
        return f"Simulation(H2O={self.n_H2O}, H3O={self.n_H3O}, Temp={self.temperature}K, Machine={self.machine})"

    def create_folder_and_molecule(self):
        """Create simulation folder and molecule structure."""
        # Versioned directory creation
        path = Path(self.folder_path)
        base = path
        i = 2
        while path.exists():
            path = Path(f"{base}_{i}")
            i += 1
        path.mkdir(parents=True)
        self.folder_path = path

        # Prepare structure and tagged xyz
        structure = Structure(self.mof_file)
        structure.pack_solvent(self.h2o_file, self.n_H2O)
        structure.pack_solvent(self.h3o_file, self.n_H3O)
        structure.name_output_file()
        structure.write_regions_and_charge()

        # Move tagged xyz into sim folder
        output_file = Path(structure.tagged_file)
        target_file = self.folder_path / output_file.name
        output_file.rename(target_file)

        # Write run.run inside the folder
        self.write_run_file_plams(output_file.name)

        log(f"Created simulation folder: {self.folder_path}")

    def write_run_file_plams(self, molecule_filename):
        """Generate run.run file using PLAMS Settings."""
        s = Settings()
        s.input.ams.Task = "MolecularDynamics"

        # MD settings
        md = s.input.ams.MolecularDynamics
        md.NSteps = self.n_steps
        md.TimeStep = 1.0
        md.Trajectory.SamplingFreq = self.sampling_freq

        # Initial velocities
        md.InitialVelocities.Type = "Random"
        md.InitialVelocities.Temperature = self.temperature

        # Thermostat
        md.Thermostat.Type = "NHC"
        md.Thermostat.Temperature = self.temperature
        md.Thermostat.Tau = 100

        # Barostat
        md.Barostat.Type = "none"
        md.Barostat.Pressure = 100000.0
        md.Barostat.Tau = 2500

        # Geometry file
        s.input.ams.System.GeometryFile = molecule_filename

        # MLPotential engine
        s.input["MLPotential"].Model = "M3GNet-UP-2022"
        s.input["MLPotential"].NumThreads = self.sim_threads

        # Prepare dummy job and extract input
        job = AMSJob(settings=s, name="temp_job")
        ams_input = job.get_input()

        run_path = self.folder_path / "run.run"

        # Write bash wrapper
        with open(run_path, "w") as f:
            f.write("#!/bin/bash\n\n")
            f.write("export NSCM={}\n\n".format(self.NSCM))
            f.write("$ADFBIN/ams << eor\n")
            f.write(ams_input)
            f.write("eor\n")

        # Make executable
        run_path.chmod(0o755)


# =====================================================================
# CLI INTERFACE
# =====================================================================
def build_parser() -> argparse.ArgumentParser:
    """Build command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Create AMS MD simulation folder with run.run file."
    )
    parser.add_argument("-N", "--num_of_mols", type=int, required=True, help="Total number of molecules")
    parser.add_argument("-H", "--H3O", type=int, required=True, help="Number of H3O molecules")
    parser.add_argument("-T", "--temperature", type=int, required=True, help="Temperature in K")
    parser.add_argument("-S", "--n_steps", type=int, default=10000, help="Number of MD steps")
    parser.add_argument("-f", "--sampling_freq", type=int, default=100, help="Sampling frequency")
    parser.add_argument("-j", "--sim_threads", type=int, default=64, help="Threads per simulation")
    parser.add_argument("--machine", type=str, default="MOPS", choices=["MOPS", "HPC", "SNELLIUS"], help="Target machine")
    parser.add_argument("--base_dir", type=str, default=os.path.join("master_thesis_project", "AMS_simulations"), help="Base directory for simulations")
    parser.add_argument("--mol_path", type=str, default=os.path.join("master_thesis_project", "AMS_simulations", "molecules"), help="Path to molecule files")
    return parser


def main(argv: Optional[list] = None) -> None:
    """Main entry point."""
    args = build_parser().parse_args(argv)

    # Setup paths
    user_home = Path.home()
    base_dir = Path(os.path.join(user_home, args.base_dir))
    mol_path = Path(os.path.join(user_home, args.mol_path))

    mof_file = mol_path / "UiO-66.xyz"
    h2o_file = mol_path / "H2O.xyz"
    h3o_file = mol_path / "H3O.xyz"

    # Validate
    for f in [mof_file, h2o_file, h3o_file]:
        if not f.exists():
            raise FileNotFoundError(f"Molecule file not found: {f}")

    # Calculate H2O count
    n_h2o = args.num_of_mols - args.H3O
    if n_h2o < 0:
        raise ValueError("H3O count cannot exceed total molecules")

    # Create simulation
    simulation = Simulation(
        base_path=base_dir,
        n_H2O=n_h2o,
        n_H3O=args.H3O,
        temperature=args.temperature,
        mof_file=mof_file,
        h2o_file=h2o_file,
        h3o_file=h3o_file,
        n_steps=args.n_steps,
        sampling_freq=args.sampling_freq,
        sim_threads=args.sim_threads,
        machine=args.machine,
    )

    simulation.create_folder_and_molecule()

    # Cleanup PLAMS directories
    for folder in Path(".").iterdir():
        if folder.is_dir() and folder.name.startswith("plams"):
            shutil.rmtree(folder, ignore_errors=True)


if __name__ == "__main__":
    main()