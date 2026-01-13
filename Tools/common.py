import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import ase.io
import ase.visualize.plot
import os
import pandas as pd
from pathlib import Path
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

# Constants
KB_EV           = 8.617333262e-5    # eV/K
KB_J            = 1.380649e-23      # J/K
H_J             = 6.62607015e-34    # J·s
H_EV            = 4.135667696e-15   # eV·s
E_CHARGE        = 1.602176634e-19   # C
ANGSTROM_TO_M   = 1e-10             # -
PS_TO_S         = 1e-12             # -


def hist2df(filename):
    return pd.read_csv(filename, sep=r"\s+").iloc[:, :3]


def rdf2df(filename):
    return pd.read_csv(filename, sep=r"\s+").iloc[:, :6]


def averagesize2df(filename, include_xxxx=False):
    """Convert averagesize.dat to DataFrame."""
    df = pd.read_csv(filename, sep=r"\s+", comment="#", 
                     names=("Group", "NumAtoms", "AverageSize", "TotalSize"))
    if include_xxxx:
        return df
    return df.drop(df[df["Group"].str.contains("xxxx", na=False)].index)


def singlegroupproperty2df(filename, prefix="col"):
    """Convert PrintProperties output to DataFrame."""
    df = pd.read_csv(filename, sep=r"\s+", header=None, skiprows=1)
    df.columns = ["#iter"] + [f"{prefix}{i+1}" for i in range(df.shape[1] - 1)]
    return df


def msd2df(filename):
    """Read MSD file to DataFrame."""
    try:
        df = pd.read_csv(filename, sep=r"\s+", skiprows=2, header=0, engine="python")
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=["#t(ps)", "value", "sum", "count"]) 

    df.columns = [str(c).strip() for c in df.columns]
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if "#t(ps)" in df.columns:
        df = df.dropna(subset=["#t(ps)"])
        df = df.sort_values(by="#t(ps)")
    return df


def tcf2df(filename):
    """Read time correlation function file to DataFrame."""
    try:
        df = pd.read_csv(filename, sep=r"\s+", skiprows=2, header=0, engine="python")
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=["#t(ps)", "value"]) 
    df.columns = [str(c).strip() for c in df.columns]
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if "#t(ps)" in df.columns:
        df = df.dropna(subset=["#t(ps)"])
        df = df.sort_values(by="#t(ps)")
    return df


def wip_seaborn_rdf_multiplot(rdf_files):
    df = pd.concat(rdf2df(filename).assign(source=filename) for filename in rdf_files)
    sns.relplot(df, x="#r", y="RDF", hue="source", kind="line", height=3)


def plot_trajectory(filename, n_images=3, plot_cell=True, rotation="-85x,-5y,0z", radii=0.6):
    """Plot a .xyz file."""
    trajectory = ase.io.read(str(filename), ":")
    n_images = min(n_images, len(trajectory))
    lattice_line_counter = -1
    if plot_cell:
        with open(filename, "r") as f:
            for line in f:
                if line.startswith("XYZ "):
                    lattice_line_counter += 1
                    splitline = line.split()
                    trajectory[lattice_line_counter].cell = np.array(
                        [float(splitline[1]), float(splitline[2]), float(splitline[3])]
                    )
    else:
        for frame in trajectory:
            frame.cell = None
    frames = np.linspace(0, len(trajectory) - 1, num=n_images, endpoint=True, dtype=np.int64)
    fig, ax = plt.subplots(1, len(frames), figsize=(10, 3))
    basename = os.path.basename(filename)
    if n_images == 1:
        ax = [ax]
    for i, frame in enumerate(frames):
        ase.visualize.plot.plot_atoms(trajectory[frame], ax=ax[i], rotation=rotation, radii=radii)
        ax[i].set_title(f"{basename} #{frame+1}")
        ax[i].axis("off")
    return ax


def donor_acceptor_cn2df(filename):
    """Parse donor_acceptor_cn file from LionAnalysis."""
    timesteps, avg_donor_cn, avg_acceptor_cn = [], [], []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            values = line.strip().split()
            if len(values) >= 3:
                timesteps.append(int(values[0]))
                data_values = [float(v) for v in values[1:]]
                mid = len(data_values) // 2
                avg_donor_cn.append(np.mean(data_values[:mid]))
                avg_acceptor_cn.append(np.mean(data_values[mid:]))
    
    return pd.DataFrame({
        'timestep': timesteps,
        'avg_donor_cn': avg_donor_cn,
        'avg_acceptor_cn': avg_acceptor_cn
    })


def double_exp(x, A, tau1, tau2):
    return A * np.exp(-x / tau1) + (1 - A) * np.exp(-x / tau2)


def single_exp(x, tau):
    return np.exp(-x / tau)


def fit_single_exp(t, y, tau_guess=1.0):
    """Fit single exponential to TCF."""
    popt, _ = curve_fit(single_exp, t, y, p0=[tau_guess], bounds=(0, np.inf))
    return popt[0], single_exp(t, popt[0])


def fit_double_exp(t, y, A_guess=0.5, tau1_guess=0.5, tau2_guess=5.0):
    """Fit double exponential to TCF."""
    popt, _ = curve_fit(
        double_exp, t, y, 
        p0=[A_guess, tau1_guess, tau2_guess],
        bounds=([0, 0, 0], [1, np.inf, np.inf])
    )
    A, tau1, tau2 = popt
    lifetime = A * tau1 + (1 - A) * tau2
    return A, tau1, tau2, lifetime, double_exp(t, A, tau1, tau2)


def extract_lifetime_from_tcf(t, y, method="single_exp"):
    """Extract lifetime τ from a time correlation function."""
    y_norm = y / y[0] if y[0] != 0 else y
    fit_info = {"method": method}
    
    if method == "single_exp":
        tau, fitted = fit_single_exp(t, y_norm)
        fit_info["tau"] = tau
        fit_info["fitted"] = fitted
    elif method == "double_exp":
        A, tau1, tau2, tau, fitted = fit_double_exp(t, y_norm)
        fit_info["A"] = A
        fit_info["tau1"] = tau1
        fit_info["tau2"] = tau2
        fit_info["lifetime"] = tau
        fit_info["fitted"] = fitted
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return tau, fit_info


def extract_ptfel_barrier(delta, free_energy, savgol_window=31, savgol_order=3, **kwargs):
    """Extract barrier from 1D PTFEL using Savitzky-Golay smoothing.
    
    Args:
        delta: Array of delta values (Å)
        free_energy: Array of free energy values (kT)
        savgol_window: Window length for SavGol filter (odd, default 31)
        savgol_order: Polynomial order for SavGol filter (default 3)
    
    Returns:
        dict with: barrier_kT, delta_saddle, fe_saddle, delta_min_left/right, 
                   fe_min_left/right, delta_min, fe_smooth, delta_sorted
    """
    from scipy.interpolate import interp1d
    
    delta = np.asarray(delta)
    fe = np.asarray(free_energy)
    
    # Sort by delta
    idx = np.argsort(delta)
    delta = delta[idx]
    fe = fe[idx]
    
    # Adjust window size
    window = min(savgol_window, len(fe) // 2 * 2 - 1)
    if window < 5:
        window = 5
    if window % 2 == 0:
        window -= 1
    
    # Apply SavGol smoothing
    fe_smooth = savgol_filter(fe, window, savgol_order)
    
    # Interpolation function for plotting
    spline_func = interp1d(delta, fe_smooth, kind='cubic', bounds_error=False, fill_value='extrapolate')
    
    # Find saddle: max in central region (-0.5 to 0.5)
    center_mask = (delta >= -0.5) & (delta <= 0.5)
    saddle_idx = np.argmax(fe_smooth[center_mask])
    delta_saddle = delta[center_mask][saddle_idx]
    fe_saddle = fe_smooth[center_mask][saddle_idx]
    
    # Find minima on left and right of saddle
    left_mask = delta < delta_saddle
    right_mask = delta > delta_saddle
    
    delta_min_left = delta[left_mask][np.argmin(fe_smooth[left_mask])]
    fe_min_left = fe_smooth[left_mask].min()
    delta_min_right = delta[right_mask][np.argmin(fe_smooth[right_mask])]
    fe_min_right = fe_smooth[right_mask].min()
    
    # Barrier = saddle height - average of minima
    fe_min_avg = (fe_min_left + fe_min_right) / 2
    barrier_kT = fe_saddle - fe_min_avg
    
    return {
        'barrier_kT': barrier_kT,
        'delta_saddle': delta_saddle,
        'fe_saddle': fe_saddle,
        'delta_min_left': delta_min_left,
        'fe_min_left': fe_min_left,
        'delta_min_right': delta_min_right,
        'fe_min_right': fe_min_right,
        'delta_min': (abs(delta_min_left) + abs(delta_min_right)) / 2,
        'spline_func': spline_func,
        'fe_smooth': fe_smooth,
        'delta_sorted': delta,
    }
