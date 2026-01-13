# Lion Analysis GUI - User Guide

A Streamlit-based graphical user interface for trajectory analysis and visualization.

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Launching the Application](#launching-the-application)
- [Features](#features)
  - [Run Analysis](#run-analysis)
  - [Interactive Plots](#interactive-plots)
- [Usage Examples](#usage-examples)
- [Troubleshooting](#troubleshooting)

## Overview

The Lion Analysis GUI provides an intuitive interface for:
- **Running trajectory analysis** with configurable parameters
- **Batch processing** multiple trajectories in parallel
- **Interactive data visualization** with customization options
- **Exporting results** in various formats (PNG, SVG, CSV)

## Installation

### Prerequisites

Ensure you have the following installed:
- Python 3.8 or higher
- AMS/SCM software suite (for `amspython`)
- Required Python packages

### Install Dependencies

```bash
# Install Streamlit and required packages
pip install streamlit plotly pandas numpy kaleido
```

**Note:** If using `amspython`, you may need to install packages in the AMS Python environment:

```bash
# Using amspython's pip
amspython -m pip install streamlit plotly pandas numpy kaleido
```

### Package Versions
- `streamlit >= 1.28.0`
- `plotly >= 5.14.0`
- `pandas >= 1.5.0`
- `numpy >= 1.23.0`
- `kaleido >= 0.2.1` (for image export)

## Launching the Application

### Using Streamlit

Navigate to the Tools directory and run:

```bash
cd Tools/
streamlit run streamlit_app.py
```

The application will open in your default web browser at `http://localhost:8501`

### Using amspython

If you need AMS integration:

```bash
cd Tools/
amspython -m streamlit run streamlit_app.py
```

### Custom Port

To use a different port:

```bash
streamlit run streamlit_app.py --server.port 8502
```

### Remote Server

To access from another machine:

```bash
streamlit run streamlit_app.py --server.address 0.0.0.0
```

Then access via `http://SERVER_IP:8501`

## Features

### 🆕 Session Persistence (NEW)

**State Management Across Pages**
- All user inputs and selections are preserved when navigating between pages
- Filter selections, analysis parameters, and plot configurations persist
- Running analyses are not interrupted when switching pages
- Analysis logs and results remain accessible

**What is Preserved:**
- File selections and filter settings
- Analysis flags and parameters (stride, threads, overwrite mode)
- Plot configurations (titles, labels, scales, colors)
- Loaded data and analysis outputs

**Benefits:**
- Seamless workflow: switch between running analysis and plotting results
- No need to reconfigure settings after page navigation
- Multi-tasking: set up multiple analyses or plots in different tabs

### Run Analysis

The Run Analysis page allows you to execute trajectory analysis with the following options:

#### 1. Analysis Modes

**Single Trajectory Mode**
- Analyze one trajectory file at a time
- Real-time log output
- Immediate feedback

**Batch Processing Mode**
- Process multiple trajectories in parallel
- Uses the `zanalyse_trajectories.sh` script
- Configurable parallel job count
- Progress tracking

#### 2. File Selection

**Browse with Filters** (NEW)
- Select from available trajectory files in `~/master_thesis_project/Trajectories/`
- **Advanced filtering** with multi-select checkboxes:
  - **Directory**: Filter by subdirectory (e.g., `xyz/`, `xyz_MOF/`)
  - **H2O Count**: Filter by number of water molecules (e.g., 29, 33, 65)
  - **H3O Count**: Filter by number of hydronium ions (e.g., 0, 4)
  - **Temperature**: Filter by temperature (e.g., 300K, 350K, 400K)
  - **Sampling Frequency**: Filter by sampling interval (e.g., 50fs, 100fs, 200fs)
- Filters combine with OR within each category and AND across categories
- **Session Persistence**: Filter selections are preserved when navigating between pages

**Enter Path Manually**
- Specify custom file path
- Useful for files outside default directories
- Path is preserved in session state

**Batch Directory**
- Enter directory path containing multiple `*_traj.xyz` files
- Shows count of discovered files

#### 3. Analysis Types

Select one or more analysis types:

| Analysis Type | Flag | Description |
|--------------|------|-------------|
| **RDF** | `--rdf` | Radial Distribution Functions |
| **MSD Water** | `--msd_water` | Mean Square Displacement (Water) |
| **MSD Hydronium** | `--msd_hydronium` | Mean Square Displacement (Hydronium) |
| **H-Bonds** | `--h_bonds` | Hydrogen Bond Analysis |
| **PTFEL** | `--ptfel` | Proton Transfer Events |
| **Residence Times** | `--res_times` | Ion-Site Residence Times |
| **Atom XYZ** | `--atom_xyz` | Track Specific Atoms |

#### 4. Parameters

**Stride**
- Frame interval for analysis
- Higher values = faster analysis, lower time resolution
- Default: 1 (analyze every frame)

**Threads (Single Mode)**
- Number of parallel threads
- Default: 8
- Increase for faster analysis on multi-core systems

**Threads per Job (Batch Mode)**
- Threads allocated to each trajectory
- Balance with max jobs for optimal performance

**Max Parallel Jobs (Batch Mode)**
- Number of trajectories processed simultaneously
- Default: 40
- Adjust based on system resources

**Overwrite Mode**
- `skip`: Skip existing results
- `prompt`: Ask before overwriting (CLI only, use `skip` or `force` in GUI)
- `force`: Always overwrite existing results

#### 5. Atom Tracking

When "Track Specific Atoms" is enabled:
- Enter atom indices (space or comma separated)
- Example: `1 2 3` or `1,2,3`
- Tracks XYZ positions over time

#### 6. Output

**Single Mode**
- **Multi-source logging panel** with tabs (NEW):
  - **📄 .out File**: Real-time display of LionAnalysis output file
  - **🐍 Python Output**: Captured stdout/stderr from analysis
  - **📁 Analysis Logs**: Browse log files from `Analysis_logs/` directory
- Success/failure status
- Download buttons for log files
- **Session Persistence**: Analysis logs remain accessible after page navigation

**Batch Mode**
- Live output streaming
- Progress updates
- Summary statistics
- Log files saved to `~/master_thesis_project/LionAnalysis/Analysis_logs/`

### Interactive Plots

The Interactive Plots page provides powerful data visualization:

#### 1. File Selection

**Browse Files**
- Scans `~/master_thesis_project/LionAnalysis/` for data files
- Filter by path (e.g., `data/`, `msd/`)
- Multi-select support

**Upload Files**
- Upload `.txt` or `.dat` files directly
- Supports multiple files

**Enter Paths**
- Manually enter file paths (one per line)
- Useful for custom locations

#### 2. Data Loading

- Automatic delimiter detection (space, tab, comma)
- Skips comment lines (starting with `#`)
- Shows row/column count for each file
- Error handling with informative messages

#### 3. Plot Types

- **Line Plot**: Continuous lines with optional markers
- **Scatter Plot**: Data points only
- **Bar Plot**: Bar chart visualization
- **Smart defaults** (NEW): Plot type is auto-suggested based on detected analysis type

#### 4. Series Customization

For each data series:

**Column Selection** (Enhanced with Smart Defaults)
- **Intelligent column detection** (NEW): Automatically selects appropriate columns based on file type
  - RDF files: Distance vs g(r)
  - MSD files: Time vs MSD
  - TCF files: Time vs Correlation
  - Histogram files: Value vs Probability
  - CN files: Time vs Coordination Number
- Manual override available for all selections

**Appearance**
- **Color**: Color picker for custom colors
- **Legend Label**: Custom label for legend
- **Line Style**: solid, dash, dot, dashdot
- **Line Width**: 1-5 pixels
- **Marker**: none, circle, square, diamond, cross, x, triangle
- **Marker Size**: 1-10 pixels

#### 5. Plot Styling

**Title and Labels** (Enhanced with Smart Defaults)
- **Auto-populated titles and labels** (NEW): Based on detected analysis type
  - RDF: "Radial Distribution Function", axes: "Distance (Å)", "g(r)"
  - MSD: "Mean Square Displacement", axes: "Time (ps)", "MSD (Å²)"
  - TCF: "Time Correlation Function", axes: "Time (ps)", "Correlation"
- **Session Persistence** (NEW): All customizations are preserved when navigating pages
- Font size control (10-24)

**Axis Configuration**
- **Scale**: Linear or logarithmic
- **Range**: Auto or manual (custom min/max)

**Legend**
- Show/hide legend
- Position: top right, top left, bottom right, bottom left

**Grid**
- Show/hide grid
- Grid style: solid, dash, dot

#### 6. Export Options

**PNG Export**
- High-resolution (300 DPI)
- 1200x800 pixels
- Publication quality

**SVG Export**
- Vector format
- Scalable without quality loss
- Ideal for manuscripts

**CSV Export**
- Combined dataset with all series
- Column headers for each series
- Easy data sharing

## Usage Examples

### Example 1: Analyze Single Trajectory

1. Open the GUI
2. Navigate to "Run Analysis"
3. Select "Single Trajectory" mode
4. Browse or enter path to `.xyz` file
5. Check "RDF" and "MSD Water"
6. Set stride to 10, threads to 8
7. Click "Run Analysis"
8. Monitor progress in log output

### Example 2: Batch Process Multiple Trajectories

1. Navigate to "Run Analysis"
2. Select "Batch Processing" mode
3. Enter directory path (e.g., `~/master_thesis_project/Trajectories/xyz/MD_130H2O/`)
4. Check desired analysis types
5. Set parameters:
   - Stride: 5
   - Threads per job: 8
   - Max parallel jobs: 40
6. Click "Run Batch Analysis"
7. Monitor overall progress

### Example 3: Plot MSD Data

1. Navigate to "Interactive Plots"
2. Select "Browse Files" mode
3. Filter by path: `msd/`
4. Select 2-3 MSD files
5. Click "Load Data"
6. For each series:
   - Set X column: time
   - Set Y column: MSD
   - Choose unique color
   - Set legend label (e.g., "300K", "350K", "400K")
7. Configure plot:
   - Title: "Mean Square Displacement"
   - X-axis: "Time (ps)"
   - Y-axis: "MSD (Å²)"
   - Enable grid
8. Click "Download PNG" to save

### Example 4: Compare RDF at Different Temperatures

1. Load multiple RDF files
2. Select "Line Plot"
3. Customize each series with different colors
4. Set X-axis label: "Distance (Å)"
5. Set Y-axis label: "g(r)"
6. Enable legend
7. Export as SVG for publication

## Troubleshooting

### Application Won't Start

**Issue**: `streamlit: command not found`

**Solution**: Install Streamlit or use full path:
```bash
python -m streamlit run streamlit_app.py
```

### Import Errors

**Issue**: `ModuleNotFoundError: No module named 'trajectory'`

**Solution**: Ensure you're running from the `Tools/` directory:
```bash
cd Tools/
streamlit run streamlit_app.py
```

### Trajectory Import Failed

**Issue**: Analysis fails with import error

**Solution**: 
- Verify `trajectory.py` is in `Tools/` directory
- Check that `analysis_inputs.py` is present
- Ensure AMS environment is activated if using `amspython`

### Batch Script Not Found

**Issue**: `zanalyse_trajectories.sh not found`

**Solution**:
- Ensure script is in `~/bin/` or `~/master_thesis_project/bin/`
- Make script executable: `chmod +x zanalyse_trajectories.sh`
- Check script path in error message

### Plot Export Failed

**Issue**: PNG/SVG export doesn't work

**Solution**: Install kaleido:
```bash
pip install kaleido
# or
amspython -m pip install kaleido
```

### Data File Loading Failed

**Issue**: "Could not parse file"

**Solution**:
- Check file format (should be space/tab/comma delimited)
- Verify file isn't empty
- Check for proper column structure
- Try opening in text editor to inspect format

### Slow Performance

**Issue**: GUI is slow or unresponsive

**Solutions**:
- Reduce number of parallel jobs in batch mode
- Increase stride for faster analysis
- Close unused browser tabs
- Check system resources (CPU, memory)

### Port Already in Use

**Issue**: `Port 8501 is already in use`

**Solution**: Use different port:
```bash
streamlit run streamlit_app.py --server.port 8502
```

## Advanced Configuration

### Streamlit Configuration

Create `~/.streamlit/config.toml`:

```toml
[server]
port = 8501
address = "localhost"
maxUploadSize = 200

[theme]
primaryColor = "#FF6B35"
backgroundColor = "#FFFFFF"
secondaryBackgroundColor = "#F0F2F6"
textColor = "#262730"
font = "sans serif"
```

### Custom File Locations

Edit `gui/utils.py` to change default paths:

```python
def get_trajectory_files(base_path=None):
    if base_path is None:
        base_path = Path("/custom/path/to/trajectories")
    # ...
```

## Tips and Best Practices

1. **Start Small**: Test with single trajectory before batch processing
2. **Monitor Resources**: Watch CPU/memory usage during batch jobs
3. **Use Stride**: For quick tests, use larger stride values
4. **Use Filters** (NEW): Leverage the advanced file filters to quickly find trajectories by temperature, composition, etc.
5. **Session Persistence** (NEW): Take advantage of preserved state—start an analysis on one page and check plots on another
6. **Smart Defaults** (NEW): Let the GUI auto-detect analysis types and populate plot settings
7. **Regular Exports**: Export plots as both PNG (viewing) and SVG (publication)
8. **Log Files**: Check log files in `Analysis_logs/` for detailed output—now accessible via tabbed interface
9. **Overwrite Wisely**: Use `skip` mode to avoid recomputing existing results

## Support

For issues or questions:
- Check the main [README.md](README.md)
- Review trajectory.py documentation
- Check AMS/SCM documentation
- Contact: [Your Email/Institution]

## Version History

- **v1.1.0** (2025-12-13): Enhanced UI with smart features
  - **Advanced file filtering** with multi-select checkboxes (H2O, H3O, temperature, sampling)
  - **Session state persistence** across page navigation
  - **Real-time .out file monitoring** with multi-source log display
  - **Smart plot defaults** with auto-detection of analysis types
  - Intelligent column selection for RDF, MSD, TCF, histogram, and CN data
  - Auto-populated axis labels and titles based on file type
  - Thread management to prevent analysis interruption during navigation

- **v1.0.0** (2025-12-13): Initial release
  - Single trajectory analysis
  - Batch processing support
  - Interactive plotting
  - Multiple export formats

## License

[Same as main project]
