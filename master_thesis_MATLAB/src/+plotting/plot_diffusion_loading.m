function fh = plot_diffusion_loading(Diffusion_data, opts)
%PLOT_DIFFUSION_LOADING Plot diffusion coefficient vs loading (H2O+H3O).
%  fh = PLOT_DIFFUSION_LOADING(Diffusion_data, Name=Value)
%
%  Produces lines of D(loading) for different temperatures. Each line
%  corresponds to a single temperature; x-axis is the "loading" defined
%  as H2O + H3O. If multiple simulations exist for the same (loading,T)
%  the mean D is plotted and a vertical errorbar shows the sample std.
%
%  Inputs:
%    Diffusion_data : 1byN struct array with fields .H2O, .H3O, .temperature, .D
%
%  Name-Value pairs (opts):
%    Scale          : "linear" (default) or "log" - applied to y-axis (D)
%    Save           : logical false (default) - save generated figure
%    SaveDir        : directory to write saved figure ("figures")
%    MarkerSize     : numeric marker size (default 6)
%    LineWidth      : numeric line width (default 1.5)
%    LegendLocation : legend placement (default "southwest")
%    FigSize        : Figure size [width, height] (default: [800, 400])
%    FontSize       : Font size for labels (default: 14)
%    TitleFontSize  : Font size for title (default: 16)
%    Verbose        : print messages (default true)
%
%  Axis behaviour added / notes:
%    - X axis (loading) is fixed to [0, 135].
%    - For linear y-scale the y-axis minimum is set to 0.
%    - For log y-scale a sensible positive ymin is chosen from data.
%
%  Output:
%    fh : figure handle (empty if nothing plotted)

    arguments
        Diffusion_data (1,:) struct
        opts.Scale (1,1) string = "log"
        opts.Save (1,1) logical = false
        opts.SaveDir (1,1) string = "figures"
        opts.MarkerSize (1,1) double = 6
        opts.LineWidth (1,1) double = 1.5
        opts.LegendLocation (1,1) string = "southwest"
        opts.FigSize (1,2) double = [800, 400]
        opts.FontSize (1,1) double = 16
        opts.TitleFontSize (1,1) double = 18
        opts.Verbose (1,1) logical = true
    end

    % Basic sanity checks
    if isempty(Diffusion_data)
        if opts.Verbose, warning('Empty Diffusion_data input. Nothing to plot.'); end
        fh = [];
        return
    end
    Table = struct2table(Diffusion_data);

    if ~ismember('D', Table.Properties.VariableNames)
        error('plot_diffusion_loading:MissingD', 'Diffusion_data must contain field .D (diffusion coefficients).');
    end
    if ~ismember('H2O', Table.Properties.VariableNames) || ~ismember('H3O', Table.Properties.VariableNames)
        error('plot_diffusion_loading:MissingFields', 'Diffusion_data must contain fields .H2O and .H3O.');
    end

    % Compute loading (H2O + H3O)
    Loading = Table.H2O + Table.H3O;

    % Keep only rows with finite D and temperature
    valid = isfinite(Table.D) & isfinite(Table.temperature);
    if ~any(valid)
        error('plot_diffusion_loading:NoValidData', 'No valid (finite) D and temperature values found.');
    end
    Table = Table(valid, :);
    Loading = Loading(valid);

    % Unique sorted loadings and temperatures
    unique_loadings = unique(Loading);
    unique_temps    = unique(Table.temperature);
    n_loadings = numel(unique_loadings);
    n_temps = numel(unique_temps);

    % Preallocate stats arrays: rows = temps, cols = loadings
    meanD = nan(n_temps, n_loadings);
    stdD  = nan(n_temps, n_loadings);
    npts  = zeros(n_temps, n_loadings);

    for ti = 1:n_temps
        temp_val = unique_temps(ti);
        mask_temp = Table.temperature == temp_val;
        for li = 1:n_loadings
            load_val = unique_loadings(li);
            % find rows matching this temp and loading
            rows = mask_temp & ( (Table.H2O + Table.H3O) == load_val );
            Ds = Table.D(rows);
            if ~isempty(Ds)
                meanD(ti,li) = mean(Ds);
                stdD(ti,li)  = std(Ds);
                npts(ti,li)  = numel(Ds);
            end
        end
    end

    % Create figure
    fh = figure('Name','Diffusion vs Loading', ...
    	'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]); 
    hold on; grid on; box on;
    ax = gca;
    col = ax.ColorOrder;
    ncol = size(col,1);

    % Plot one line per temperature. Errorbars when npts>1
    for ti = 1:n_temps
        y = meanD(ti, :);
        x = unique_loadings;

        % Only plot points where meanD is finite
        mask = isfinite(y);
        if ~any(mask)
            if opts.Verbose
                fprintf('No data for %g K — skipping\n', unique_temps(ti));
            end
            continue
        end

        c = col(mod(ti-1, ncol) + 1, :);

        % Plot line + markers for available points
        plot(x(mask), y(mask), '-o', ...
            'Color', c, ...
            'MarkerFaceColor', c, ...
            'MarkerSize', opts.MarkerSize, ...
            'LineWidth', opts.LineWidth, ...
            'DisplayName', sprintf('%g K', unique_temps(ti)));

        % Add errorbars (std) if multiple points at that (loading, temp)
        err = stdD(ti, :);
        % Only draw errorbars where std is finite and > 0 and corresponding mean is finite
        errmask = isfinite(err) & (err > 0) & mask;
        if any(errmask)
            % errorbar draws separate markers by default; suppress marker by LineStyle 'none'
            errorbar(x(errmask), y(errmask), err(errmask), 'LineStyle','none', ...
                    'Color', c, 'HandleVisibility','off', ...
                    'LineWidth', opts.LineWidth * 0.8);
        end
    end

    % Labels and title
    xlabel('Loading (H2O/Cell)', 'FontSize', opts.FontSize);
    ylabel('D [m^2 s^{-1}]', 'FontSize', opts.FontSize);
    title('Diffusion coefficient vs loading', ... 
          'FontSize', opts.TitleFontSize, 'FontWeight', 'bold');
    ax.FontSize = opts.FontSize;

    % Set Y-scale
    try
        set(ax, 'YScale', opts.Scale);
    catch
        warning('Invalid Scale value. Using ''linear''.');
        set(ax, 'YScale', 'linear');
        opts.Scale = "linear";
    end

    % Axis limits
    axis padded
    leg = legend('Location', opts.LegendLocation);
    leg.FontSize = opts.FontSize - 1;

    % Save if requested
    if opts.Save
        if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
        out = fullfile(opts.SaveDir, 'diffusion_vs_loading.png');
        exportgraphics(ax, out, 'Resolution', 300);
        if opts.Verbose, fprintf('Saved: %s\n', out); end
    end

    % Return only if figure is still valid
    if ~isgraphics(fh)
        fh = [];
    end
end