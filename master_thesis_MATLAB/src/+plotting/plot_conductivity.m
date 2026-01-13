function fh = plot_conductivity(Conductivity_data, opts)
%PLOT_CONDUCTIVITY Plot proton conductivity vs temperature and loading.
%  fh = PLOT_CONDUCTIVITY(Conductivity_data, Name=Value)
%
%  Creates plots showing vehicular, structural, and total conductivity
%  with consistent loading colors matching RDF plots.
%
%  Inputs:  
%    Conductivity_data : table with Loading, Temperature, sigma_veh, sigma_struct, sigma_total
%
%  Name-Value options:  
%    PlotType      : "vs_temp" | "vs_loading" | "both" (default: "vs_temp")
%    Scale         : "linear" | "log" (default:  "log")
%    Save          : Save figures (default: false)
%    SaveDir       : Save directory (default: "figures/conductivity")
%    Verbose       : Print messages (default: true)
%    
%    FigSize       : Figure size [width, height] (default: [800, 400])
%    FontSize      : Font size (default: 14)
%    TitleFontSize : Title font size (default: 16)
%    LineWidth     : Line width (default:  1.5)
%    MarkerSize    : Marker size (default: 6)

    arguments
        Conductivity_data (:,11) table
        
        opts.PlotType (1,1) string = "vs_temp"
        opts.Scale (1,1) string = "log"
        opts.Save (1,1) logical = false
        opts.SaveDir (1,1) string = "figures/conductivity"
        opts.Verbose (1,1) logical = true
        
        opts.FigSize (1,2) double = [800, 400]
        opts.FontSize (1,1) double = 14
        opts.TitleFontSize (1,1) double = 16
        opts.LineWidth (1,1) double = 1.5
        opts.MarkerSize (1,1) double = 6
    end
    
    % Get unique values
    unique_temps = unique(Conductivity_data.Temperature);
    unique_loadings = unique(Conductivity_data.Loading);
    
    fh = [];
    
    % Plot vs temperature
    if ismember(opts.PlotType, ["vs_temp", "both"])
        fh_temp = plot_vs_temperature(Conductivity_data, unique_temps, unique_loadings, opts);
        fh = [fh; fh_temp];
    end
    
    % Plot vs loading
    if ismember(opts.PlotType, ["vs_loading", "both"])
        fh_load = plot_vs_loading(Conductivity_data, unique_temps, unique_loadings, opts);
        fh = [fh; fh_load];
    end
end

%% Helper: Plot Arrhenius (ln(sigma*T) vs 1/T)
function fh = plot_vs_temperature(Table, unique_temps, unique_loadings, opts)
    % Initialize figure
    fh = figure('Name', 'Arrhenius Plot: ln(sigma*T) vs 1/T', ...
                'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]);
    ax = gca; hold on; grid on; box on;
    
    % Color Logic
    ax_temp = axes('Visible', 'off');
    default_colors = ax_temp.ColorOrder;
    delete(ax_temp); 
    loading_colors = [default_colors(4:end,:); [1.0, 0.6, 0.0]];
    
    h_loading = []; 
    
    for i = 1:numel(unique_loadings)
        loading = unique_loadings(i);
        mask = Table.Loading == loading;
        
        % Transformations
        T = Table.Temperature(mask);
        invT_1000 = 1000 ./ T;
        
        c = loading_colors(mod(i-1, size(loading_colors, 1)) + 1, :);
        
        % Helper to calculate ln(sigma * 0.01 * T)
        calc_lnST = @(sigma) log((sigma .* 0.01) .* T);
        
        % Plot Vehicular
        y_veh = calc_lnST(Table.sigma_veh(mask));
        plot(invT_1000, y_veh, ':', 'Color', c, 'LineWidth', opts.LineWidth, ...
             'MarkerSize', opts.MarkerSize, 'HandleVisibility', 'off');
        
        % Plot Structural
        if ismember('sigma_struct', Table.Properties.VariableNames)
            y_struct = calc_lnST(Table.sigma_struct(mask));
            if any(isfinite(y_struct))
                plot(invT_1000, y_struct, '--', 'Color', c, 'LineWidth', opts.LineWidth, ...
                     'MarkerSize', opts.MarkerSize, 'HandleVisibility', 'off');
            end
        end
        
        % Plot Total
        y_total = calc_lnST(Table.sigma_total(mask));
        h_loading(i) = plot(invT_1000, y_total, '-', 'Color', c, 'LineWidth', opts.LineWidth, ...
             'MarkerSize', opts.MarkerSize, 'DisplayName', sprintf('%d H_2O', loading));
    end
    
    % Mechanism Legend
    h_veh = plot(nan, nan, 'k:', 'DisplayName', 'Vehicular');
    h_str = plot(nan, nan, 'k--', 'DisplayName', 'Structural');
    h_tot = plot(nan, nan, 'k-', 'DisplayName', 'Total');
    
    legend([h_loading, h_veh, h_str, h_tot], 'Location', 'northeastoutside', 'FontSize', opts.FontSize - 2);
    
    % Labels and Formatting
    xlabel('1000/T [K^{-1}]', 'FontSize', opts.FontSize);
    ylabel('ln(\sigma T / [S K cm^{-1}])', 'FontSize', opts.FontSize);
    title('Arrhenius Analysis of Proton Conductivity', 'FontSize', opts.TitleFontSize, 'FontWeight', 'bold');
    
    set(ax, 'YScale', 'linear', 'FontSize', opts.FontSize);
    
    if opts.Save
        if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
        exportgraphics(ax, fullfile(opts.SaveDir, 'arrhenius_log_sigmaT_cm.png'), 'Resolution', 300);
    end
end

%% Helper: Plot vs loading
function fh = plot_vs_loading(Table, unique_temps, unique_loadings, opts)
    % FIX: Fresh figure
    fh = figure('Name', 'Conductivity vs Loading', ...
                'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]);
    ax = gca; hold on; grid on; box on;
    
    temp_colors = [0, 0, 1; 0, 0.6, 0; 1, 0, 0]; % Blue, Green, Red
    h_temp = [];
    
    for i = 1:numel(unique_temps)
        temp = unique_temps(i);
        mask = Table.Temperature == temp;
        c = temp_colors(min(i, size(temp_colors, 1)), :);
        
        % Vehicular (Solid)
        plot(Table.Loading(mask), Table.sigma_veh(mask), 'o-', 'Color', c, ...
             'LineWidth', opts.LineWidth, 'MarkerSize', opts.MarkerSize, 'HandleVisibility', 'off');
        
        % Total (Dashed) - Represents the Temp in legend
        h_temp(i) = plot(Table.Loading(mask), Table.sigma_total(mask), 's--', 'Color', c, ...
             'LineWidth', opts.LineWidth, 'MarkerSize', opts.MarkerSize, 'DisplayName', sprintf('%d K', temp));
    end
    
    % Dummy entries for mechanism styles
    h_m1 = plot(nan, nan, 'ko-', 'DisplayName', 'Vehicular');
    h_m2 = plot(nan, nan, 'ks--', 'DisplayName', 'Total');
    
    legend([h_temp, h_m1, h_m2], 'Location', 'northeastoutside', 'FontSize', opts.FontSize - 2);
    
    xlabel('Loading [H_2O/Cell]', 'FontSize', opts.FontSize);
    ylabel('Conductivity [S/m]', 'FontSize', opts.FontSize);
    title('Proton Conductivity vs Loading', 'FontSize', opts.TitleFontSize, 'FontWeight', 'bold');
    set(ax, 'YScale', opts.Scale, 'FontSize', opts.FontSize);
    
    if opts.Save
        if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
        exportgraphics(ax, fullfile(opts.SaveDir, 'conductivity_vs_loading.png'), 'Resolution', 300);
    end
end