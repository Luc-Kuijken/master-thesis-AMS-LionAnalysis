function opts = style_applier(preset, plot_type, base_opts)
%STYLE_APPLIER Apply preset styles for all plotting functions.
%  opts = STYLE_APPLIER(preset, plot_type) returns Name-Value options.
%  opts = STYLE_APPLIER(preset, plot_type, base_opts) merges with existing options.
%
%  Presets:
%    "individual"  : Full featured for interactive viewing
%    "subfigure"   : Minimal - for LaTeX subfigures (no title, no legend)
%    "publication" : Optimized for journal figures
%    "presentation": Larger fonts for slides
%
%  Plot types:
%    "msd"         : MSD vs time plots
%    "arrhenius"   : ln(D) vs 1/T plots
%    "beta"        : Beta exponent vs time plots
%    "diffusion"   : Diffusion vs loading plots
%    "heatmap"     : 2D/3D heatmap plots
%    "general"     : Generic settings for any plot
%
%  Example:
%    opts = plotting.style_applier("publication", "msd");
%    plotting.plot_msd(Diffusion_data, Mol_values, opts{:});

    arguments
        preset (1,1) string = "individual"
        plot_type (1,1) string = "general"
        base_opts (1,:) cell = {}
    end

    % Validate inputs
    valid_presets = ["individual", "subfigure", "publication", "presentation"];
    valid_types = ["msd", "arrhenius", "beta", "diffusion", "heatmap", "general"];
    
    if ~ismember(lower(preset), valid_presets)
        error('Unknown preset: "%s".Use: %s', preset, strjoin(valid_presets, ', '));
    end
    if ~ismember(lower(plot_type), valid_types)
        error('Unknown plot_type: "%s". Use: %s', plot_type, strjoin(valid_types, ', '));
    end

    % Get base style from preset
    base_style = get_preset_style(preset);
    
    % Get plot-specific adjustments
    plot_style = get_plot_style(plot_type, preset);
    
    % Merge: plot_style overrides base_style
    style = merge_structs(base_style, plot_style);
    
    % Convert struct to cell array of Name-Value pairs
    opts = struct_to_cell(style);
    
    % Merge with base options (base_opts override everything)
    if ~isempty(base_opts)
        opts = [opts, base_opts];
    end
end

%% Helper: Get base style from preset
function style = get_preset_style(preset)
    switch lower(preset)
        case "individual"
            style = struct(...
                'ShowTitle', true, ...
                'ShowLegend', true, ...
                'ShowColorbar', true, ...
                'FontSize', 12, ...
                'LabelFontSize', 14, ...
                'TitleFontSize', 14, ...
                'LegendFontSize', 11, ...
                'AxisLineWidth', 1.0, ...
                'FigureWidth', 560, ...
                'FigureHeight', 420, ...
                'DPI', 300, ...
                'DisplayMode', "interactive", ...
                'Verbose', true ...
            );
            
        case "subfigure"
            style = struct(...
                'ShowTitle', false, ...
                'ShowLegend', false, ...
                'ShowColorbar', true, ...
                'FontSize', 14, ...
                'LabelFontSize', 16, ...
                'TitleFontSize', 16, ...
                'LegendFontSize', 14, ...
                'AxisLineWidth', 0.8, ...
                'FigureWidth', 520, ...
                'FigureHeight', 420, ...
                'AxisWidth', 380, ...
                'AxisHeight', 380, ...
                'DPI', 300, ...
                'DisplayMode', "off", ...
                'Verbose', true ...
            );
            
        case "publication"
            style = struct(...
                'ShowTitle', false, ...
                'ShowLegend', true, ...
                'ShowColorbar', true, ...
                'FontSize', 10, ...
                'LabelFontSize', 11, ...
                'TitleFontSize', 12, ...
                'LegendFontSize', 9, ...
                'AxisLineWidth', 0.8, ...
                'FigureWidth', 560, ...
                'FigureHeight', 480, ...
                'AxisWidth', 400, ...
                'AxisHeight', 400, ...
                'DPI', 600, ...
                'DisplayMode', "interactive", ...
                'Verbose', true ...
            );
            
        case "presentation"
            style = struct(...
                'ShowTitle', true, ...
                'ShowLegend', true, ...
                'ShowColorbar', true, ...
                'FontSize', 16, ...
                'LabelFontSize', 18, ...
                'TitleFontSize', 20, ...
                'LegendFontSize', 16, ...
                'AxisLineWidth', 1.2, ...
                'FigureWidth', 700, ...
                'FigureHeight', 600, ...
                'DPI', 600, ...
                'DisplayMode', "interactive", ...
                'Verbose', true ...
            );
    end
end

%% Helper: Get plot-specific style adjustments
function style = get_plot_style(plot_type, preset)
    style = struct();
    
    switch lower(plot_type)
        case "msd"
            style.LegendLocation = "northwest";
            if preset == "publication"
                style.FigureHeight = 300;
            end
            
        case "arrhenius"
            style.LegendLocation = "northeast";
            if preset == "publication"
                style.MarkerSize = 5;
            elseif preset == "presentation"
                style.MarkerSize = 10;
            end
            
        case "beta"
            style.LegendLocation = "southeast";
            style.LineAlpha = 0.2;
            style.AvgLineWidth = 1.5;
            if preset == "presentation"
                style.AvgLineWidth = 2.5;
            end
            
        case "diffusion"
            style.LegendLocation = "southwest";
            style.Scale = "log";
            
        case "heatmap"
            % Heatmap-specific defaults
            style.HeatmapAlpha = 0.4;
            style.AtomAlpha = 0.4;
            style.MarkerSize = 3;
            style.MaxValue = 15;
            style.AxisWidth = 400;
            style.AxisHeight = 400;
            style.FigureWidth = 560;
            style.FigureHeight = 480;
            
        case "general"
            style.LegendLocation = "best";
    end
end

%% Helper: Merge two structs (second overrides first)
function merged = merge_structs(base, override)
    merged = base;
    if isempty(override)
        return
    end
    fields = fieldnames(override);
    for i = 1:numel(fields)
        merged.(fields{i}) = override.(fields{i});
    end
end

%% Helper: Convert struct to cell array of Name-Value pairs
function opts = struct_to_cell(style)
    fields = fieldnames(style);
    opts = cell(1, 2*numel(fields));
    for i = 1:numel(fields)
        opts{2*i-1} = fields{i};
        opts{2*i} = style.(fields{i});
    end
end