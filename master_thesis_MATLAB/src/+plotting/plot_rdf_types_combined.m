function fh = plot_rdf_types_combined(RDF_data, rdf_types_to_plot, opts)
%PLOT_RDF_TYPES_COMBINED Plot specific RDF types for all loadings in 2by2 layout.
%  fh = PLOT_RDF_TYPES_COMBINED(RDF_data, rdf_types_to_plot, Name=Value)
%
%  Creates one figure with 2by2 subplots, one per loading (33, 65, 98, 130 H2O).
%  Each subplot shows lines for the specified RDF types at different temperatures.
%
%  Inputs:
%    RDF_data          : struct array from io.load_rdf_struct
%    rdf_types_to_plot : string array of RDF types to plot
%                        e.g., ["LINKEROH", "LINKEROO"] or ["MU3OH", "MU3OO"]
%
%  Name-Value options:
%    XLim              : x-axis limits (default: [1, 7])
%    YLim              : y-axis limits (default: auto)
%    ShowPMF           : Show PMF instead of RDF (default: false)
%    Save              : Save figure (default: false)
%    SaveDir           : Save directory (default: "figures/rdf")
%    LineWidth         : Line width (default:  1.0)
%    Verbose           : Print messages (default: true)
%    
%    FigSize           : Figure size [width, height] (default: [1000, 600])
%    FontSize          : Font size (default: 12)
%    TitleFontSize     : Title font size (default: 14)
%    
%    HorizontalSpacing : Spacing between columns (default: 0.02)
%    VerticalSpacing   : Spacing between rows (default: 0.04)
%    MarginLeft        : Left margin (default: 0.08)
%    MarginRight       : Right margin (default:  0.02)
%    MarginBottom      : Bottom margin (default: 0.10)
%    MarginTop         : Top margin (default: 0.08)
%    
%    LineStyles        : Cell array of line styles for each RDF type
%                        (default: {'-', '--', '-.', ':'})
%    FigureTitle       : Override default title (default: auto-generated)
%    Filename          : Override default filename (default:  auto-generated)

    arguments
        RDF_data (1,: ) struct
        rdf_types_to_plot (:,1) string
        
        opts.XLim (1,2) double = [1, 7]
        opts.YLim (1,: ) double = []
        opts.ShowPMF (1,1) logical = false
        opts.Save (1,1) logical = false
        opts.SaveDir (1,1) string = "figures/rdf"
        opts.LineWidth (1,1) double = 1.0
        opts.Verbose (1,1) logical = true
        
        opts.FigSize (1,2) double = [1000, 600]
        opts.FontSize (1,1) double = 12
        opts.TitleFontSize (1,1) double = 14
        
        opts.HorizontalSpacing (1,1) double = 0.02
        opts.VerticalSpacing (1,1) double = 0.04
        opts.MarginLeft (1,1) double = 0.08
        opts.MarginRight (1,1) double = 0.02
        opts.MarginBottom (1,1) double = 0.10
        opts.MarginTop (1,1) double = 0.08
        
        opts.LineStyles (1,:) cell = {'-', '--', '-.', ': '}
        opts.FigureTitle (1,1) string = ""
        opts.Filename (1,1) string = ""
    end
    
    % Convert to table for easier filtering
    Table = struct2table(RDF_data);
    Table.Total = Table.H2O + Table.H3O;
    
    % Filter to only requested RDF types
    type_mask = ismember(Table.rdf_type, rdf_types_to_plot);
    if ~any(type_mask)
        error('No data found for RDF types: %s', strjoin(rdf_types_to_plot, ', '));
    end
    Table = Table(type_mask, :);
    RDF_data_filtered = RDF_data(type_mask);
    
    % Get unique loadings and temperatures
    unique_loadings = sort(unique(Table.Total));
    unique_temps = sort(unique(Table.temperature));
    
    n_loadings = numel(unique_loadings);
    n_types = numel(rdf_types_to_plot);
    
    % Define colors for temperatures
    temp_colors = [
        0, 0,   1;      % Blue for 300K
        0, 0.6, 0;      % Green for 350K
        1, 0,   0       % Red for 400K
    ];
    
    % Auto-generate title if not provided
    if opts.FigureTitle == ""
        if opts.ShowPMF
            opts.FigureTitle = sprintf('Potential of Mean Force - %s', ...
                                      strjoin(rdf_types_to_plot, ' & '));
        else
            opts.FigureTitle = sprintf('Radial Distribution Function - %s', ...
                                      strjoin(rdf_types_to_plot, ' & '));
        end
    end
    
    % Auto-generate filename if not provided
    if opts.Filename == ""
        type_str = lower(strjoin(rdf_types_to_plot, '_'));
        if opts.ShowPMF
            opts.Filename = sprintf('pmf_%s_all_loadings.png', type_str);
        else
            opts.Filename = sprintf('rdf_%s_all_loadings.png', type_str);
        end
    end
    
    % Set y-axis label
    if opts.ShowPMF
        y_label = 'PMF [k_BT]';
    else
        y_label = 'g(r)';
    end
    
    % Create figure
    fh = figure('Name', opts.FigureTitle, ...
                'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]);
    
    % Calculate subplot dimensions
    plot_width = (1 - opts.MarginLeft - opts.MarginRight - opts.HorizontalSpacing) / 2;
    plot_height = (1 - opts.MarginBottom - opts.MarginTop - opts.VerticalSpacing) / 2;
    
    % Define positions for each subplot
    positions = [
        opts.MarginLeft, opts.MarginBottom + plot_height + opts.VerticalSpacing, plot_width, plot_height;
        opts.MarginLeft + plot_width + opts.HorizontalSpacing, opts.MarginBottom + plot_height + opts.VerticalSpacing, plot_width, plot_height;
        opts.MarginLeft, opts.MarginBottom, plot_width, plot_height;
        opts.MarginLeft + plot_width + opts.HorizontalSpacing, opts.MarginBottom, plot_width, plot_height
    ];
    
    % Store axes handles
    ax_handles = gobjects(min(n_loadings, 4), 1);
    
    % Create subplots
    for i = 1:min(n_loadings, 4)
        loading = unique_loadings(i);
        
        ax = axes('Position', positions(i, : ));
        ax_handles(i) = ax;
        hold on; grid on; box on;
        
        % Add loading label
        text(0.02, 0.95, sprintf('%d H_2O', loading), ...
             'Units', 'normalized', ...
             'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'top', ...
             'FontSize', opts.FontSize, ...
             'FontWeight', 'bold', ...
             'BackgroundColor', 'white', ...
             'EdgeColor', [0.7 0.7 0.7], ...
             'Margin', 2);
        
        % Plot each combination of rdf_type and temperature
        for ti = 1:numel(unique_temps)
            temp = unique_temps(ti);
            c = temp_colors(min(ti, size(temp_colors, 1)), :);
            
            for ri = 1:n_types
                rdf_type = rdf_types_to_plot(ri);
                
                % Select line style based on RDF type
                line_style = opts.LineStyles{mod(ri-1, numel(opts.LineStyles)) + 1};
                
                % Find all matching simulations
                mask = (Table.Total == loading) & ...
                       (Table.temperature == temp) & ...
                       (Table.rdf_type == rdf_type);
                
                matching_idx = find(mask);
                
                if isempty(matching_idx)
                    continue
                end
                
                % Collect matching RDF structs
                matching = RDF_data_filtered(matching_idx);
                
                % Average RDF across different H3O counts
                [r_avg, y_avg] = average_rdf(matching, opts.ShowPMF);
                
                % Only show legend for subplot 2 (top right)
                if i == 2
                    display_name = sprintf('%s %dK', rdf_type, temp);
                    handle_vis = 'on';
                else
                    display_name = '';
                    handle_vis = 'off';
                end
                
                % Plot averaged curve
                plot(r_avg, y_avg, ...
                     'Color', c, ...
                     'LineStyle', line_style, ...
                     'LineWidth', opts.LineWidth, ...
                     'DisplayName', display_name, ...
                     'HandleVisibility', handle_vis);
                
                if opts.Verbose && i == 1
                    fprintf('Plotting:  loading=%d, %dK, %s (averaged %d simulations)\n', ...
                            loading, temp, rdf_type, numel(matching));
                end
            end
        end
        
        % Set axis limits
        xlim(opts.XLim);
        if ~isempty(opts.YLim)
            ylim(opts.YLim);
        end
        
        % Add reference line
        if ~opts.ShowPMF
            yline(1, ':k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        else
            yline(0, ':k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        
        % Apply font size
        ax.FontSize = opts.FontSize;
        
        % Only show legend on subplot 2 (top right)
        if i == 2
            legend('Location', 'east', 'FontSize', opts.FontSize - 2);
        end
    end
    
    % Link axes
    linkaxes(ax_handles, 'x');
    if n_loadings >= 2
        linkaxes(ax_handles(1:2), 'y');
    end
    if n_loadings >= 4
        linkaxes(ax_handles(3:4), 'y');
    end
    
    % Remove tick labels from inner axes
    if n_loadings >= 1
        ax_handles(1).XTickLabel = [];
    end
    if n_loadings >= 2
        ax_handles(2).XTickLabel = [];
        ax_handles(2).YTickLabel = [];
    end
    if n_loadings >= 4
        ax_handles(4).YTickLabel = [];
    end
    
    % Add shared x-axis label
    annotation('textbox', [0, 0, 1, opts.MarginBottom * 0.5], ...
               'String', 'r [Å]', ...
               'FontSize', opts.FontSize + 2, ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'middle', ...
               'FitBoxToText', 'off');
    
    % Add shared y-axis label
    ax_ylabel = axes('Position', [0, 0, opts.MarginLeft * 0.5, 1], 'Visible', 'off');
    text(0.5, 0.5, y_label, ...
         'Parent', ax_ylabel, ...
         'Units', 'normalized', ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', ...
         'FontSize', opts.FontSize + 2, ...
         'Rotation', 90);
    
    % Add overall title
    annotation('textbox', [0, 1 - opts.MarginTop * 0.9, 1, opts.MarginTop * 0.9], ...
               'String', opts.FigureTitle, ...
               'FontSize', opts.TitleFontSize + 2, ...
               'FontWeight', 'bold', ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'middle', ...
               'FitBoxToText', 'off');
    
    % Save if requested
    if opts.Save
        if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
        filepath = fullfile(opts.SaveDir, opts.Filename);
        exportgraphics(fh, filepath, 'Resolution', 300);
        if opts.Verbose
            fprintf('Saved: %s\n', filepath);
        end
    end
end

%% Helper:  Average RDF across multiple simulations
function [r_avg, y_avg] = average_rdf(rdf_structs, use_pmf)
    r_ref = rdf_structs(1).data.r;
    n_points = numel(r_ref);
    n_sims = numel(rdf_structs);
    
    y_all = NaN(n_points, n_sims);
    
    for k = 1:n_sims
        data = rdf_structs(k).data;
        r_k = data.r;
        if use_pmf
            y_k = data.pmf;
        else
            y_k = data.rdf;
        end
        
        if isequal(size(r_k), size(r_ref)) && all(abs(r_k - r_ref) < 1e-6)
            y_all(:, k) = y_k;
        else
            y_all(:, k) = interp1(r_k, y_k, r_ref, 'linear', NaN);
        end
    end
    
    r_avg = r_ref;
    y_avg = mean(y_all, 2, 'omitnan');
end