function fh = plot_ptfel_1d(PTFEL_data, opts)
%PLOT_PTFEL_1D Plot 1D PT-FEL with water (blue) and MOF (red) overlaid.
%  fh = PLOT_PTFEL_1D(PTFEL_data, Name=Value)
%
%  Creates one figure per (loading, temperature) combination.
%  Shows water (blue) and MOF (red) PT-FEL on the same plot.
%  Raw data plotted with alpha=0.4, smoothed data with alpha=1.0.
%
%  Inputs:
%    PTFEL_data :   struct array from io.load_ptfel_struct
%
%  Name-Value options: 
%    XLim          :  x-axis limits (default:  [-1.5, 1.5])
%    YLim          : y-axis limits (default: [0, 12])
%    Save          : Save figure (default: false)
%    SaveDir       : Save directory (default: "figures/ptfel")
%    Verbose       : Print messages (default: true)
%    
%    FigSize       : Figure size [width, height] (default: [800, 500])
%    FontSize      : Font size (default: 12)
%    TitleFontSize : Title font size (default: 14)
%    LineWidth     : Line width for smoothed (default: 2.5)
%    RawAlpha      :  Transparency for raw data (default: 0.4)

    arguments
        PTFEL_data (1,: ) struct
        
        opts.XLim (1,2) double = [-1.5, 1.5]
        opts.YLim (1,2) double = [0, 12]
        opts.Save (1,1) logical = false
        opts.SaveDir (1,1) string = "figures/ptfel"
        opts.Verbose (1,1) logical = true
        
        opts.FigSize (1,2) double = [800, 500]
        opts.FontSize (1,1) double = 12
        opts.TitleFontSize (1,1) double = 14
        opts.LineWidth (1,1) double = 2.5
        opts.RawAlpha (1,1) double = 0.4
    end
    
    % Filter for 1D data only
    Table = struct2table(PTFEL_data);
    mask_1d = ~[PTFEL_data.is_2d]';
    PTFEL_1d = PTFEL_data(mask_1d);
    Table_1d = Table(mask_1d, :);
    
    if isempty(PTFEL_1d)
        warning('No 1D PT-FEL data found.');
        fh = [];
        return
    end
    
    % Get unique (loading, temperature) combinations
    unique_combos = unique(Table_1d(: , {'H2O', 'H3O', 'temperature'}), 'rows');
    n_combos = height(unique_combos);
    
    % Colors:  blue for water, red for MOF
    color_water = [0, 0, 1];
    color_mof = [1, 0, 0];
    
    % Preallocate figure handles
    fh = gobjects(n_combos, 1);
    
    % Loop through each (loading, temperature) combination
    for i = 1:n_combos
        H2O_val = unique_combos.H2O(i);
        H3O_val = unique_combos.H3O(i);
        T_val = unique_combos.temperature(i);
        loading = H2O_val + H3O_val;
        
        if opts.Verbose
            fprintf('Plotting:  %dH2O_%dH3O, %dK\n', H2O_val, H3O_val, T_val);
        end
        
        % Create figure
        fig_name = sprintf('1D PTFEL - %d H2O - %d K', loading, T_val);
        fh(i) = figure('Name', fig_name, ...
                      'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]);
        ax = gca; hold on; grid on; box on;
        
        % Find water PTFEL for this loading and temperature
        idx_water = find(([PTFEL_1d.type] == "water") & ...
                        ([PTFEL_1d.H2O] == H2O_val) & ...
                        ([PTFEL_1d.H3O] == H3O_val) & ...
                        ([PTFEL_1d.temperature] == T_val));
        
        % Find MOF PTFEL for this loading and temperature
        idx_mof = find(([PTFEL_1d.type] == "mof") & ...
                      ([PTFEL_1d.H2O] == H2O_val) & ...
                      ([PTFEL_1d.H3O] == H3O_val) & ...
                      ([PTFEL_1d.temperature] == T_val));
        
        % Plot water PTFEL if available
        if ~isempty(idx_water)
            data_water = PTFEL_1d(idx_water);
            delta_w = data_water.data.delta;
            fe_w = data_water.data.free_energy;
            fe_w = fe_w - min(fe_w);  % Shift to zero minimum
            
            % Plot raw water data (alpha 0.4)
            plot(delta_w, fe_w, '-', 'Color', [color_water, opts.RawAlpha], ...
                 'LineWidth', 1, 'HandleVisibility', 'off');
            
            % Plot smoothed water data (alpha 1.0) if available
            if isfield(data_water, 'fe_smooth') && ~isempty(data_water.fe_smooth)
                fe_smooth = data_water.fe_smooth - min(data_water.fe_smooth);
                plot(delta_w, fe_smooth, '-', 'Color', color_water, ...
                     'LineWidth', opts.LineWidth, ...
                     'DisplayName', 'Water-Water');
            end
        end
        
        % Plot MOF PTFEL if available
        if ~isempty(idx_mof)
            data_mof = PTFEL_1d(idx_mof);
            delta_m = data_mof.data.delta;
            fe_m = data_mof.data.free_energy;
            fe_m = fe_m - min(fe_m);  % Shift to zero minimum
            
            % Plot raw MOF data (alpha 0.4)
            plot(delta_m, fe_m, '-', 'Color', [color_mof, opts.RawAlpha], ...
                 'LineWidth', 1, 'HandleVisibility', 'off');
            
            % Plot smoothed MOF data (alpha 1.0) if available
            if isfield(data_mof, 'fe_smooth') && ~isempty(data_mof.fe_smooth)
                fe_smooth = data_mof.fe_smooth - min(data_mof.fe_smooth);
                plot(delta_m, fe_smooth, '-', 'Color', color_mof, ...
                     'LineWidth', opts.LineWidth, ...
                     'DisplayName', 'MOF-Water');
            end
        end
        
        % Labels and formatting
        xlabel('\delta (Å)', 'FontSize', opts.FontSize, 'FontWeight', 'bold');
        ylabel('Free Energy (kT)', 'FontSize', opts.FontSize, 'FontWeight', 'bold');
        title(fig_name, 'FontSize', opts.TitleFontSize, 'FontWeight', 'bold');
        
        xlim(opts.XLim);
        ylim(opts.YLim);
        
        ax.FontSize = opts.FontSize;
        legend('Location', 'northwest', 'FontSize', opts.FontSize);
        
        % Save if requested
        if opts.Save
            if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
            filename = sprintf('ptfel_1d_%dH2O_%dH3O_%dK.png', H2O_val, H3O_val, T_val);
            filepath = fullfile(opts.SaveDir, filename);
            exportgraphics(ax, filepath, 'Resolution', 300);
            if opts.Verbose, fprintf('Saved: %s\n', filepath); end
        end
    end
end