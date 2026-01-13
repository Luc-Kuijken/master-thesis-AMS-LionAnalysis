function fh = plot_ptfel_1d(PTFEL_data, opts)
%PLOT_PTFEL_1D Plot 1D PT-FEL with water and MOF overlaid. 
%  fh = PLOT_PTFEL_1D(PTFEL_data, Name=Value)
%
%  Shows water (blue) and MOF (red) PT-FEL on the same plot.
%  Raw data plotted with alpha=0.4, smoothed data with alpha=1.0.
%
%  Inputs:
%    PTFEL_data :  struct array from io. load_ptfel_struct
%
%  Name-Value options:
%    Temperature   : Filter by temperature (default:  350)
%    XLim          : x-axis limits (default: [-1.5, 1.5])
%    YLim          : y-axis limits (default:  [0, 12])
%    Save          : Save figure (default: false)
%    SaveDir       : Save directory (default: "figures/ptfel")
%    Verbose       : Print messages (default: true)
%    
%    FigSize       : Figure size [width, height] (default: [1200, 400])
%    FontSize      : Font size (default: 14)
%    TitleFontSize :  Title font size (default: 16)
%    LineWidth     : Line width for smoothed (default: 2. 5)

    arguments
        PTFEL_data (1,: ) struct
        
        opts. Temperature (1,1) double = 350
        opts.XLim (1,2) double = [-1.5, 1.5]
        opts.YLim (1,2) double = [0, 12]
        opts.Save (1,1) logical = false
        opts.SaveDir (1,1) string = "figures/ptfel"
        opts.Verbose (1,1) logical = true
        
        opts.FigSize (1,2) double = [1200, 400]
        opts.FontSize (1,1) double = 14
        opts.TitleFontSize (1,1) double = 16
        opts.LineWidth (1,1) double = 2.5
    end
    
    % Filter for 1D data at specified temperature
    Table = struct2table(PTFEL_data);
    mask = ~[PTFEL_data.is_2d]' & (Table.temperature == opts. Temperature);
    PTFEL_1d = PTFEL_data(mask);
    
    if isempty(PTFEL_1d)
        warning('No 1D PT-FEL data found for T=%dK', opts.Temperature);
        fh = [];
        return
    end
    
    % Get unique loadings
    loadings = unique(Table.H2O(mask) + Table.H3O(mask));
    
    % Colors:  blue for water, red for MOF
    color_water = [0, 0, 1];
    color_mof = [1, 0, 0];
    
    % Create figure
    fig_name = sprintf('1D PTFEL - T = %d K', opts.Temperature);
    fh = figure('Name', fig_name, ...
               'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]);
    ax = gca; hold on; grid on; box on;
    
    % Plot each loading
    for i = 1:numel(loadings)
        loading = loadings(i);
        
        % Get water PT-FEL
        idx_water = find(([PTFEL_1d.type] == "water") & ...
                        (([PTFEL_1d. H2O] + [PTFEL_1d. H3O]) == loading));
        
        if ~isempty(idx_water)
            data_water = PTFEL_1d(idx_water);
            delta_w = data_water.data.delta;
            fe_w = data_water.data.free_energy - min(data_water.data.free_energy);
            
            % Plot raw water data (alpha 0.4)
            p1 = plot(delta_w, fe_w, '-', 'Color', [color_water, 0.4], ... 
                     'LineWidth', 1, 'HandleVisibility', 'off');
            
            % Plot smoothed water data (alpha 1.0)
            if isfield(data_water, 'spline_func') && ~isempty(data_water. spline_func)
                delta_fine = linspace(min(delta_w), max(delta_w), 300);
                fe_smooth = data_water.spline_func(delta_fine);
                fe_smooth = fe_smooth - min(fe_smooth);
                
                plot(delta_fine, fe_smooth, '-', 'Color', color_water, ...
                     'LineWidth', opts.LineWidth, ... 
                     'DisplayName', sprintf('%d H_2O/Cell (\\Delta G^{\\ddagger}=%.1f kT)', ...
                                           loading, data_water.barrier_kT));
            end
        end
        
        % Get MOF PT-FEL
        idx_mof = find(([PTFEL_1d. type] == "mof") & ...
                      (([PTFEL_1d.H2O] + [PTFEL_1d.H3O]) == loading));
        
        if ~isempty(idx_mof)
            data_mof = PTFEL_1d(idx_mof);
            delta_m = data_mof.data. delta;
            fe_m = data_mof.data.free_energy - min(data_mof.data.free_energy);
            
            % Plot raw MOF data (alpha 0.4)
            plot(delta_m, fe_m, '-', 'Color', [color_mof, 0.4], ...
                 'LineWidth', 1, 'HandleVisibility', 'off');
            
            % Plot smoothed MOF data (alpha 1.0)
            if isfield(data_mof, 'spline_func') && ~isempty(data_mof.spline_func)
                delta_fine = linspace(min(delta_m), max(delta_m), 300);
                fe_smooth = data_mof.spline_func(delta_fine);
                fe_smooth = fe_smooth - min(fe_smooth);
                
                % Calculate average barrier for MOF (left + right) / 2
                if isfield(data_mof, 'barrier_kT')
                    avg_barrier = data_mof.barrier_kT;
                else
                    avg_barrier = NaN;
                end
                
                plot(delta_fine, fe_smooth, '-', 'Color', color_mof, ... 
                     'LineWidth', opts.LineWidth, ...
                     'DisplayName', sprintf('%d H_2O/Cell MOF (\\Delta G^{\\ddagger}=%.1f kT)', ...
                                           loading, avg_barrier));
            end
        end
    end
    
    % Labels and formatting
    xlabel('\delta (Å)', 'FontSize', opts.FontSize, 'FontWeight', 'bold');
    ylabel('Free Energy (kT)', 'FontSize', opts.FontSize, 'FontWeight', 'bold');
    title(fig_name, 'FontSize', opts.TitleFontSize, 'FontWeight', 'bold');
    
    xlim(opts.XLim);
    ylim(opts.YLim);
    
    ax.FontSize = opts.FontSize;
    legend('Location', 'northwest', 'FontSize', opts.FontSize - 2);
    
    % Add legend annotation for line styles
    text(0.98, 0.02, 'Water-Water | MOF-Water', ... 
         'Units', 'normalized', ... 
         'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'bottom', ...
         'FontSize', opts.FontSize - 2, ...
         'Color', [0.3, 0.3, 0.3], ...
         'FontWeight', 'bold');
    
    % Save if requested
    if opts.Save
        if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
        filename = sprintf('ptfel_1d_%dK.png', opts.Temperature);
        filepath = fullfile(opts.SaveDir, filename);
        exportgraphics(ax, filepath, 'Resolution', 300);
        if opts.Verbose, fprintf('Saved: %s\n', filepath); end
    end
end