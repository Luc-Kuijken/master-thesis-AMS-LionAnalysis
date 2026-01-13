function fh = plot_rdf_combined(RDF_data, rdf_types_to_plot, opts)
%PLOT_RDF_COMBINED Plot RDF curves in a flexible layout.
%  fh = PLOT_RDF_COMBINED(RDF_data, rdf_types_to_plot, Name=Value)
%
%  Two modes:  
%    GroupBy="loading" (default): 2by2 layout, one subplot per loading, lines colored by temperature
%    GroupBy="temperature": 1by3 layout, one subplot per temperature, lines colored by loading
%
%  Each subplot can show multiple RDF types (e.g., ["O-H", "O-O"] or ["LinkerO-H", "LinkerO-O"])
%  with different line styles.
%
%  Inputs:
%    RDF_data           : struct array from io.load_rdf_struct
%    rdf_types_to_plot  : string array of RDF types for display
%                         e.g., ["O-H", "O-O"] or ["LinkerO-H", "LinkerO-O"] or ["Mu3O-H", "Mu3O-O"]
%
%  Name-Value options:
%    GroupBy           : "loading" | "temperature" (default:   "loading")
%    XLim              : x-axis limits (default: [1, 7])
%    YLim              : y-axis limits (default: auto)
%    ShowPMF           : Show PMF instead of RDF (default: false)
%    Save              : Save figure (default: false)
%    SaveDir           : Save directory (default: "figures/rdf")
%    LineWidth         : Line width (default: 1.0)
%    Verbose           : Print messages (default: true)
%    
%    FigSize           : Figure size [width, height] (default: auto-sized based on GroupBy)
%    FontSize          : Font size (default:  12)
%    TitleFontSize     : Title font size (default: 14)
%    
%    HorizontalSpacing : Spacing between columns (default: 0.02)
%    VerticalSpacing   : Spacing between rows (default: 0.04)
%    MarginLeft        : Left margin (default: 0.08)
%    MarginRight       : Right margin (default: 0.02)
%    MarginBottom      : Bottom margin (default: 0.10)
%    MarginTop         : Top margin (default: 0.08)
%    
%    LineStyles        : Cell array of line styles for RDF types (default: {'-', '--'})
%    FigureTitle       : Override default title (default: auto-generated)
%    Filename          : Override default filename (default: auto-generated)

    arguments
        RDF_data (1,:) struct
        rdf_types_to_plot (:,1) string
        
        opts.GroupBy (1,1) string = "loading"
        opts.XLim (1,2) double = [1, 7]
        opts.YLim (1,:) double = []
        opts.ShowPMF (1,1) logical = false
        opts.Save (1,1) logical = false
        opts.SaveDir (1,1) string = "figures/rdf"
        opts.LineWidth (1,1) double = 1.0
        opts.Verbose (1,1) logical = true
        
        opts.FigSize (1,2) double = [0, 0]  % Auto-size if [0,0]
        opts.FontSize (1,1) double = 12
        opts.TitleFontSize (1,1) double = 14
        
        opts.HorizontalSpacing (1,1) double = 0.02
        opts.VerticalSpacing (1,1) double = 0.04
        opts.MarginLeft (1,1) double = 0.08
        opts.MarginRight (1,1) double = 0.02
        opts.MarginBottom (1,1) double = 0.10
        opts.MarginTop (1,1) double = 0.08
        
        opts.LineStyles (1,:) cell = {'-', '--'}
        opts.FigureTitle (1,1) string = ""
        opts.Filename (1,1) string = ""
    end
    
    % Validate GroupBy
    if ~ismember(opts.GroupBy, ["loading", "temperature"])
        error('GroupBy must be "loading" or "temperature"');
    end
    
    % Convert to table for easier filtering
    Table = struct2table(RDF_data);
    Table.Total = Table.H2O + Table.H3O;
    
    % Filter to only requested RDF types
    type_mask = ismember(Table.rdf_type, rdf_types_to_plot);
    if ~any(type_mask)
        error('No data found for RDF types: %s\nAvailable types: %s', ...
              strjoin(rdf_types_to_plot, ', '), ...
              strjoin(unique(Table.rdf_type), ', '));
    end
    Table = Table(type_mask, :);
    RDF_data_filtered = RDF_data(type_mask);
    
    % Get unique values for grouping
    unique_loadings = sort(unique(Table.Total));
    unique_temps = sort(unique(Table.temperature));
    n_types = numel(rdf_types_to_plot);
    
    % Define consistent colors
    % Temperature colors (for loading-grouped plots)
    temp_colors = [
        0, 0,   1;      % Blue for 300K
        0, 0.6, 0;      % Green for 350K
        1, 0,   0       % Red for 400K
    ];
    
    % Loading colors (for temperature-grouped plots)
    ax_temp = axes('Visible', 'off');
    default_colors = ax_temp.ColorOrder;
    close(gcf);
    
    loading_colors = [default_colors(4:  end,:  ); [1.0, 0.6, 0.0]];
    
    % Auto-generate title and filename
    type_str = strjoin(rdf_types_to_plot, ' & ');
    
    if opts.FigureTitle == ""
        if opts.GroupBy == "loading"
            base_title = sprintf('RDF - %s (by Loading)', type_str);
        else
            base_title = sprintf('RDF - %s (by Temperature)', type_str);
        end
        opts.FigureTitle = base_title;
    end
    
    if opts.Filename == ""
        type_file = lower(strjoin(strrep(rdf_types_to_plot, '-', ''), '_'));
        if opts.GroupBy == "loading"
            opts.Filename = sprintf('rdf_%s_by_loading.png', type_file);
        else
            opts.Filename = sprintf('rdf_%s_by_temperature.png', type_file);
        end
    end
    
    % Set y-axis label
    if opts.ShowPMF
        y_label = 'PMF [k_BT]';
    else
        y_label = 'g(r)';
    end
    
    % Detect if this is a MOF-water interaction plot (needs YLim constraint)
    rdf_label = get_rdf_type_label(rdf_types_to_plot);
    is_mof_water = contains(lower(rdf_label), ["linkero", "mu3o"]);
    
    % Determine subplot organization and layout
    if opts.GroupBy == "loading"
        subplot_values = unique_loadings;
        n_subplots = numel(unique_loadings);
        color_values = unique_temps;
        color_map = temp_colors;
        value_label = 'H_2O';
        
        % 2×2 layout for loading
        n_rows = 2;
        n_cols = 2;
        if all(opts.FigSize == 0)
            opts.FigSize = [1000, 600];
        end
    else  % GroupBy == "temperature"
        subplot_values = unique_temps;
        n_subplots = numel(unique_temps);
        color_values = unique_loadings;
        color_map = loading_colors;
        value_label = 'K';
        
        % 1×3 horizontal layout for temperature
        n_rows = 1;
        n_cols = 3;
        if all(opts.FigSize == 0)
            opts.FigSize = [1400, 500];
        end
        
        % Increase font sizes for temperature-grouped plots by 2
        opts.FontSize = opts.FontSize + 2;
        opts.TitleFontSize = opts.TitleFontSize + 2;
    end
    
    % Create figure
    fh = figure('Name', opts.FigureTitle, ...
                'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]);
    
    % Calculate subplot dimensions
    plot_width = (1 - opts.MarginLeft - opts.MarginRight - (n_cols-1)*opts.HorizontalSpacing) / n_cols;
    plot_height = (1 - opts.MarginBottom - opts.MarginTop - (n_rows-1)*opts.VerticalSpacing) / n_rows;
    
    % Generate positions for subplots
    positions = zeros(n_rows * n_cols, 4);
    idx = 1;
    for row = 1:n_rows
        for col = 1:n_cols
            left = opts.MarginLeft + (col-1)*(plot_width + opts.HorizontalSpacing);
            bottom = opts.MarginBottom + (n_rows - row)*(plot_height + opts.VerticalSpacing);
            positions(idx, : ) = [left, bottom, plot_width, plot_height];
            idx = idx + 1;
        end
    end
    
    % Store axes handles
    ax_handles = gobjects(min(n_subplots, n_rows*n_cols), 1);
    
    % Create subplots
    for i = 1:min(n_subplots, n_rows*n_cols)
        subplot_val = subplot_values(i);
        
        ax = axes('Position', positions(i, :));
        ax_handles(i) = ax;
        hold on; grid on; box on;
        
        % Add subplot label with line style info for MOF-water plots
        if opts.GroupBy == "loading"
            if is_mof_water
                label_text = sprintf('%d H_2O\n(solid=%s, dash=%s)', ...
                                   subplot_val, rdf_types_to_plot(1), rdf_types_to_plot(2));
            else
                label_text = sprintf('%d H_2O', subplot_val);
            end
        else
            if is_mof_water
                label_text = sprintf('%d K\n(solid=%s, dash=%s)', ...
                                   subplot_val, rdf_types_to_plot(1), rdf_types_to_plot(2));
            else
                label_text = sprintf('%d K', subplot_val);
            end
        end
        
        text(0.02, 0.95, label_text, ...
             'Units', 'normalized', ...
             'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'top', ...
             'FontSize', opts.FontSize - 1, ...
             'FontWeight', 'bold', ...
             'BackgroundColor', 'white', ...
             'EdgeColor', [0.7 0.7 0.7], ...
             'Margin', 2);
        
        % Plot each combination of color_value and rdf_type
        for ci = 1:numel(color_values)
            color_val = color_values(ci);
            
            % Select color using your specified indexing
            c = color_map(mod(ci-1, size(color_map, 1)) + 1, :);
            
            for ri = 1:n_types
                rdf_type = rdf_types_to_plot(ri);
                
                % Select line style
                line_style = opts.LineStyles{mod(ri-1, numel(opts.LineStyles)) + 1};
                
                % Build filter mask
                if opts.GroupBy == "loading"
                    mask = (Table.Total == subplot_val) & ...
                           (Table.temperature == color_val) & ...
                           (Table.rdf_type == rdf_type);
                    legend_label = sprintf('%dK', color_val);
                else
                    mask = (Table.temperature == subplot_val) & ...
                           (Table.Total == color_val) & ...
                           (Table.rdf_type == rdf_type);
                    legend_label = sprintf('%dH2O', color_val);
                end
                
                matching_idx = find(mask);
                
                if isempty(matching_idx)
                    continue
                end
                
                % Collect and average RDF
                matching = RDF_data_filtered(matching_idx);
                [r_avg, y_avg] = average_rdf(matching, opts.ShowPMF);
                
                % Legend placement:   top-right for loading (2×2), rightmost for temperature (1×3)
                % Only show legend for first RDF type (ri == 1) to avoid duplicates
                if ((opts.GroupBy == "loading" && i == 2) || ...
                   (opts.GroupBy == "temperature" && i == n_cols)) && ri == 1
                    display_name = legend_label;
                    handle_vis = 'on';
                else
                    display_name = '';
                    handle_vis = 'off';
                end
                
                % Plot
                plot(r_avg, y_avg, ...
                     'Color', c, ...
                     'LineStyle', line_style, ...
                     'LineWidth', opts.LineWidth, ...
                     'DisplayName', display_name, ...
                     'HandleVisibility', handle_vis);
                
                if opts.Verbose && i == 1
                    if opts.GroupBy == "loading"
                        fprintf('Plotting:   %d H2O, %dK, %s (averaged %d)\n', ...
                                subplot_val, color_val, rdf_type, numel(matching));
                    else
                        fprintf('Plotting:  %dK, %d H2O, %s (averaged %d)\n', ...
                                subplot_val, color_val, rdf_type, numel(matching));
                    end
                end
            end
        end
        
        % Set axis limits - Apply XLim first, then auto YLim, then constrain if needed
        xlim(opts.XLim);
        
        % Auto-scale Y if not specified, or apply user YLim
        if ~isempty(opts.YLim)
            ylim(opts.YLim);
        else
            % Let MATLAB auto-scale first
            drawnow;
            current_ylim = ylim;
            
            % For MOF-water interactions, limit to [0, 2] max
            if is_mof_water
                ylim([0, min(2, current_ylim(2))]);
            end
        end
        
        % Add reference line
        if ~opts.ShowPMF
            yline(1, ':  k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        else
            yline(0, ': k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        
        % Apply font size
        ax.FontSize = opts.FontSize;
        
        % Legend placement - simplified for MOF-water
        if (opts.GroupBy == "loading" && i == 2) || ...
           (opts.GroupBy == "temperature" && i == n_cols)
            
            if is_mof_water
                % For MOF-water:   compact legend showing only colors
                leg = legend('Location', 'northeast');
                leg.FontSize = opts.FontSize - 2;
                leg.Box = 'on';
                leg.Color = [1 1 1 0.95];  % Semi-transparent white
            else
                % For water-water:  original full legend
                leg = legend('Location', 'east', 'FontSize', opts.FontSize - 2);
            end
        end
    end
    
    % Link axes
    linkaxes(ax_handles, 'x');
    
    % For 2×2 layout:   link rows separately
    if opts.GroupBy == "loading"
        if n_subplots >= 2
            linkaxes(ax_handles(1:2), 'y');
        end
        if n_subplots >= 4
            linkaxes(ax_handles(3:4), 'y');
        end
        
        % Remove tick labels from inner axes
        if n_subplots >= 1, ax_handles(1).XTickLabel = []; end
        if n_subplots >= 2
            ax_handles(2).XTickLabel = [];
            ax_handles(2).YTickLabel = [];
        end
        if n_subplots >= 4, ax_handles(4).YTickLabel = []; end
    else
        % For 1×3 layout:  link all y-axes
        linkaxes(ax_handles, 'y');
        
        % Remove y-tick labels from middle and right panels
        for i = 2:min(n_subplots, n_cols)
            ax_handles(i).YTickLabel = [];
        end
    end
    
    % Shared labels
    annotation('textbox', [0, 0, 1, opts.MarginBottom * 0.5], ...
               'String', 'r [Å]', ...
               'FontSize', opts.FontSize + 2, ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'middle', ...
               'FitBoxToText', 'off');
    
    ax_ylabel = axes('Position', [0, 0, opts.MarginLeft * 0.5, 1], 'Visible', 'off');
    text(0.5, 0.5, y_label, ...
         'Parent', ax_ylabel, ...
         'Units', 'normalized', ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', ...
         'FontSize', opts.FontSize + 2, ...
         'Rotation', 90);
    
    annotation('textbox', [0, 1 - opts.MarginTop * 0.9, 1, opts.MarginTop * 0.9], ...
               'String', opts.FigureTitle, ...
               'FontSize', opts.TitleFontSize + 2, ...
               'FontWeight', 'bold', ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'middle', ...
               'FitBoxToText', 'off');
    
    % Save
    if opts.Save
        if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
        filepath = fullfile(opts.SaveDir, opts.Filename);
        exportgraphics(fh, filepath, 'Resolution', 300);
        if opts.Verbose
            fprintf('Saved:  %s\n', filepath);
        end
    end
end

%% Helper:   Get RDF type label for title
function label = get_rdf_type_label(rdf_types)
    % Extract the atom type from RDF type names
    % ["O-H", "O-O"] → "O"
    % ["LinkerO-H", "LinkerO-O"] → "LinkerO"
    % ["Mu3O-H", "Mu3O-O"] → "Mu3O"
    
    % Get first type and extract prefix before hyphen
    first_type = char(rdf_types(1));
    dash_idx = strfind(first_type, '-');
    
    if isempty(dash_idx)
        label = first_type;
    else
        label = first_type(1:dash_idx(1)-1);
    end
end

%% Helper: Average RDF
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