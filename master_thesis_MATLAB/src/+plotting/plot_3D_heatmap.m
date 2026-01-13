function fh = plot_3D_heatmap(Heatmap_data_3D, opts)
%PLOT_3D_HEATMAP Plot 3D heatmap as isosurface and/or scatter plot.
%  fh = PLOT_3D_HEATMAP(Heatmap_data_3D, Name=Value)
%  Inputs: Heatmap_data_3D struct with .data3D containing X,Y,Z,V
%  Name-Value options:
%    Index         : Which heatmap to plot (default: 1)
%    PlotType      : "isosurface" | "scatter" | "both" (default: "isosurface")
%    
%    IsoValue      : Isosurface threshold (0-1, default: 0.5)
%    FaceColor     : Isosurface color (default: 'cyan')
%    FaceAlpha     : Isosurface transparency (default: 0.7)
%    ShowSlices    : Show 2D slices with isosurface (default: false)
%    
%    MaxValue      : Max value for scatter colormap (default: max from data)
%    HeatmapAlpha  : Power for alpha mapping (default: 0.4)
%    MarkerSize    : Scatter marker size (default: 3)
%    Downsample    : Downsample factor for scatter (1=all points, default: 1)
%    
%    Colormap      : Colormap name (default: "turbo")
%    ShowColorbar  : Show colorbar (default: true)
%    
%    MOFFile       : Path to MOF .xyz file (default: "data/processed/UiO-66.xyz")
%    ShowMOF       : Show MOF unit cell (default: true)
%    AtomAlpha     : MOF atom opacity, 0=no atoms (default: 0)
%    
%    Save          : Save figure (default: false)
%    SaveDir       : Save directory (default: "figures/heatmaps")
%    Verbose       : Print messages (default: true)

    arguments
        Heatmap_data_3D (1,:) struct
        
        opts.Index (1,1) double = 1
        opts.PlotType (1,1) string = "isosurface"
        
        opts.IsoValue (1,1) double = 0.5
        opts.FaceColor = 'cyan'
        opts.FaceAlpha (1,1) double = 0.7
        opts.ShowSlices (1,1) logical = false
        
        opts.MaxValue (1,1) double = NaN
        opts.HeatmapAlpha (1,1) double = 0.4
        opts.MarkerSize (1,1) double = 3
        opts.Downsample (1,1) double = 1
        
        opts.Colormap (1,1) string = "turbo"
        opts.ShowColorbar (1,1) logical = true
        
        opts.MOFFile (1,1) string = "data/processed/UiO-66.xyz"
        opts.ShowMOF (1,1) logical = true
        opts.AtomAlpha (1,1) double = 0
        
        opts.Save (1,1) logical = false
        opts.SaveDir (1,1) string = "figures/heatmaps"
        opts.Verbose (1,1) logical = true
    end
    
    % Validate inputs
    if opts.Index < 1 || opts.Index > numel(Heatmap_data_3D)
        error('Index out of range (1 to %d)', numel(Heatmap_data_3D));
    end
    
    if ~ismember(opts.PlotType, ["isosurface", "scatter", "both"])
        error('PlotType must be "isosurface", "scatter", or "both"');
    end
    
    % Load MOF structure if requested
    if opts.ShowMOF && isfile(opts.MOFFile)
        mof = io.load_mof_structure(opts.MOFFile, Verbose=false);
    else
        mof = [];
        if opts.ShowMOF && opts.Verbose
            warning('MOF file not found: %s', opts.MOFFile);
        end
    end
    
    % Extract data
    heatmap = Heatmap_data_3D(opts.Index);
    X = heatmap.data3D.X;
    Y = heatmap.data3D.Y;
    Z = heatmap.data3D.Z;
    V = heatmap.data3D.V;
    
    % Calculate time and unit for filename
    time = heatmap.sampling_freq * 500;
    if time < 1e5
        time = time/1e3; unit = 'ps';
    else
        time = time/1e6; unit = 'ns';
    end
    
    % Determine value limits
    if isnan(opts.MaxValue)
        v_max = max(V(:), [], 'omitnan');
    else
        v_max = opts.MaxValue;
    end
    v_min = 0;
    
    % Generate figure name
    if opts.PlotType == "both"
        plot_suffix = "Combined";
    else
        plot_suffix = char(opts.PlotType);
        plot_suffix(1) = upper(plot_suffix(1));
    end
    fig_name = sprintf('3D %s - %dH2O %dH3O %dK %d%s', ...
                      plot_suffix, heatmap.H2O, heatmap.H3O, ...
                      heatmap.temperature, time, unit);
    
    % Create figure
    fh = figure('Name', fig_name);
    ax = gca; hold on; grid on; box on;
    
    % Plot based on type
    if opts.PlotType == "isosurface" || opts.PlotType == "both"
        plot_isosurface_view(X, Y, Z, V, opts, v_min, v_max, ax);
    end
    
    if opts.PlotType == "scatter" || opts.PlotType == "both"
        plot_scatter_view(X, Y, Z, V, opts, v_min, v_max, ax);
    end
    
    % Overlay MOF unit cell
    if opts.ShowMOF && ~isempty(mof)
        plotting.plot_mof_frame(mof, ...
                               Projection="3D", ...
                               LineStyle="--", ...
                               AtomScale=1.5, ...
                               AtomAlpha=opts.AtomAlpha);
    end
    
    % Colormap and colorbar (only for scatter)
    if opts.PlotType == "scatter" || opts.PlotType == "both"
        colormap(ax, opts.Colormap);
        caxis([v_min, v_max]);
        
        if opts.ShowColorbar
            cb = colorbar;
            cb.Label.String = 'Density';
        end
    end
    
    % Labels and view
    xlabel('Position X [Å]');
    ylabel('Position Y [Å]');
    zlabel('Position Z [Å]');
    title(fig_name);
    axis equal tight;
    view(135, 35);
    rotate3d on;
    
    % Save if requested
    if opts.Save
        if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
        filename = sprintf('AvOPs_%dH2O_%dH3O_%dK_%d%s_3D_%s.png', ...
                          heatmap.H2O, heatmap.H3O, heatmap.temperature, ...
                          time, unit, lower(plot_suffix));
        filepath = fullfile(opts.SaveDir, filename);
        exportgraphics(ax, filepath, 'Resolution', 300);
        if opts.Verbose, fprintf('Saved: %s\n', filepath); end
    end
end

%% Helper: Plot isosurface view
function plot_isosurface_view(X, Y, Z, V, opts, v_min, v_max, ax)
    % Normalize values to [0, 1] for isosurface
    V_norm = (V - v_min) / (v_max - v_min);
    V_norm(isnan(V_norm)) = 0;
    
    % Create isosurface
    isosurf = isosurface(X, Y, Z, V_norm, opts.IsoValue);
    
    if ~isempty(isosurf.vertices)
        p = patch(isosurf, 'Parent', ax);
        p.FaceColor = opts.FaceColor;
        p.EdgeColor = 'none';
        p.FaceAlpha = opts.FaceAlpha;
        
        % Add lighting for better 3D perception
        lighting gouraud;
        camlight('headlight');
        material dull;
    else
        warning('No isosurface found at IsoValue=%.2f', opts.IsoValue);
    end
    
    % Optionally show 2D slices
    if opts.ShowSlices
        [ny, nx, nz] = size(V);
        
        % XY slice at mid-Z
        slice(X, Y, Z, V, [], [], Z(1,1,round(nz/2)));
        
        % YZ slice at mid-X
        slice(X, Y, Z, V, X(1,round(nx/2),1), [], []);
        
        % ZX slice at mid-Y
        slice(X, Y, Z, V, [], Y(round(ny/2),1,1), []);
        
        shading interp;
        colormap(ax, opts.Colormap);
    end
end

%% Helper: Plot scatter view
function plot_scatter_view(X, Y, Z, V, opts, v_min, v_max, ax)
    % Flatten 3D grid to 1D arrays
    x = X(:);
    y = Y(:);
    z = Z(:);
    v = V(:);
    
    % Remove NaN values
    v = ((v./3).^2)-1;
    
    valid = isfinite(v) & (v > v_min);
    x = x(valid);
    y = y(valid);
    z = z(valid);
    v = v(valid);
    v(v > v_max-5) = v_max-5;
    
    % Randomized downsampling (for performance)
    if opts.Downsample >= 10
        % Compute number of samples to keep (1/N of total)
        keepCount = ceil(numel(x) / opts.Downsample);

        % Generate random permutation of indices
        randIdx = randperm(numel(x), keepCount);

        % Apply random selection
        x = x(randIdx);
        y = y(randIdx);
        z = z(randIdx);
        v = v(randIdx);
    else
        error("DownSample too small, must be 10 or above")
    end
    
    if isempty(x)
        warning('No valid data points to plot in scatter view');
        return
    end
    
    % Create scatter plot with color
    s = scatter3(x, y, z, v.^1.3, v, 'filled');
    
    % Apply alpha mapping if requested
    if opts.HeatmapAlpha < 1
        factor_3d = 0.5;  % correction factor for 3D view
        alpha_vals = calc_alpha(v, v_min, v_max, opts.HeatmapAlpha * factor_3d);
        s.AlphaData = alpha_vals;
        s.MarkerFaceAlpha = 'flat';
    end
end

%% Helper: Calculate alpha values
function alpha = calc_alpha(v, v_min, v_max, alpha_value)
    % Normalize to [0, 1] and apply power transform
    alpha = (v - v_min) / (v_max - v_min);
    alpha = max(0, min(1, alpha));
    power = 1 - alpha_value;  % inverted to make sure low alpha_value gives transparent points
    alpha = alpha.^power;
end