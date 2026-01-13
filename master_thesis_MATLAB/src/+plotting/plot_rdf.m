function fh = plot_rdf(RDF_data, rdf_types_to_plot, opts)
%PLOT_RDF Plot individual RDF figures for each loading or temperature.  
%  fh = PLOT_RDF(RDF_data, rdf_types_to_plot, Name=Value)
%
%  Creates separate figures for each unique loading or temperature value.
%  Two modes:     
%    GroupBy="loading": One figure per loading, lines show different temperatures
%    GroupBy="temperature": One figure per temperature, lines show different loadings
%
%  Inputs:
%    RDF_data          : struct array from io.load_rdf_struct
%    rdf_types_to_plot : string array of RDF types, e.g., ["O-H", "O-O"]
%
%  Name-Value options:
%    GroupBy          : "loading" | "temperature" (default: "loading")
%    XLim             : x-axis limits (default: [1, 7])
%    YLim             : y-axis limits (default: auto)
%    ShowPMF          : Show PMF instead of RDF (default: false)
%    Save             : Save figures (default: false)
%    SaveDir          : Save directory (default: "figures/rdf")
%    LineWidth        : Line width (default: 1.5)
%    LegendLocation   : Legend placement (default: "best")
%    Verbose          : Print messages (default: true)
%    FigSize          : Figure size [width, height] (default: [800, 400])
%    FontSize         : Font size for labels and ticks (default: 14)
%    TitleFontSize    : Font size for title (default: 16)

    arguments
        RDF_data (1,:) struct
        rdf_types_to_plot (:,1) string
        
        opts.GroupBy (1,1) string = "loading"
        opts.XLim (1,2) double = [1, 7]
        opts.YLim (1,:  ) double = []
        opts.ShowPMF (1,1) logical = false
        opts.Save (1,1) logical = false
        opts.SaveDir (1,1) string = "figures/rdf"
        opts.LineWidth (1,1) double = 1.5
        opts.LegendLocation (1,1) string = "best"
        opts.Verbose (1,1) logical = true
        opts.FigSize (1,2) double = [800, 400]
        opts.FontSize (1,1) double = 16
        opts.TitleFontSize (1,1) double = 18
    end
    
    % Validate GroupBy
    if ~ismember(opts.GroupBy, ["loading", "temperature"])
        error('GroupBy must be "loading" or "temperature"');
    end
    
    % Convert to table
    Table = struct2table(RDF_data);
    Table.Total = Table.H2O + Table.H3O;
    
    % Filter to requested RDF types
    type_mask = ismember(Table.rdf_type, rdf_types_to_plot);
    if ~any(type_mask)
        error('No data found for RDF types: %s', strjoin(rdf_types_to_plot, ', '));
    end
    Table = Table(type_mask, :);
    RDF_data_filtered = RDF_data(type_mask);
    
    % Get unique values
    unique_loadings = sort(unique(Table.Total));
    unique_temps = sort(unique(Table.temperature));
    n_types = numel(rdf_types_to_plot);
    
    % Define colors
    temp_colors = [
        0, 0,   1;      % Blue for 300K
        0, 0.6, 0;      % Green for 350K
        1, 0,   0       % Red for 400K
    ];
    
    % Loading colors
    ax_temp = axes('Visible', 'off');
    default_colors = ax_temp.ColorOrder;
    close(gcf);
    loading_colors = [default_colors(4:end,:); [1.0, 0.6, 0.0]];
    
    % Set y-axis label
    if opts.ShowPMF
        y_label = 'PMF [k_BT]';
        plot_type = 'PMF';
    else
        y_label = 'g(r)';
        plot_type = 'RDF';
    end
    
    % Generate RDF type label for title (remove "-", e.g., "O-H" becomes "O")
    rdf_type_label = get_rdf_type_label(rdf_types_to_plot);
    
    % Detect if this is a MOF-water interaction plot (needs YLim constraint)
    is_mof_water = any(contains(lower(rdf_types_to_plot), ["linkero", "mu3o"]));
    
    % Determine grouping
    if opts.GroupBy == "loading"
        group_values = unique_loadings;
        color_values = unique_temps;
        color_map = temp_colors;
    else
        group_values = unique_temps;
        color_values = unique_loadings;
        color_map = loading_colors;
    end
    
    n_groups = numel(group_values);
    fh = gobjects(n_groups, 1);
    
    % Create one figure per group
    for i = 1:n_groups
        group_val = group_values(i);
        
        % Generate figure name with RDF type
        if opts.GroupBy == "loading"
            fig_name = sprintf('%s - %s - %d H2O', plot_type, rdf_type_label, group_val);
            filename_base = sprintf('%s_%s_%dH2O', lower(plot_type), lower(strrep(rdf_type_label, '-', '')), group_val);
        else
            fig_name = sprintf('%s - %s - %d K', plot_type, rdf_type_label, group_val);
            filename_base = sprintf('%s_%s_%dK', lower(plot_type), lower(strrep(rdf_type_label, '-', '')), group_val);
        end
        
        % Create figure with custom size
        fh(i) = figure('Name', fig_name, ...
                      'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]);
        ax = gca; hold on; grid on; box on;
        
        xlabel('r [Å]', 'FontSize', opts.FontSize);
        ylabel(y_label, 'FontSize', opts.FontSize);
        title_suffix = ' (solid=O-H, dashed=O-O)';
        title(sprintf('%s%s', fig_name, title_suffix), 'FontSize', opts.TitleFontSize, 'FontWeight', 'bold');
        ax.FontSize = opts.FontSize;
        
        % Plot each combination of color_value and rdf_type
        for ci = 1:numel(color_values)
            color_val = color_values(ci);
            c = color_map(mod(ci-1, size(color_map, 1)) + 1, :);
            
            for ri = 1:n_types
                rdf_type = rdf_types_to_plot(ri);
                
                % Line style: solid for first type, dashed for second
                if ri == 1
                    line_style = '-';
                else
                    line_style = '--';
                end
                
                % Build filter mask
                if opts.GroupBy == "loading"
                    mask = (Table.Total == group_val) & ...
                           (Table.temperature == color_val) & ...
                           (Table. rdf_type == rdf_type);
                    if ri == 1
                        legend_label = sprintf('%dK', color_val);
                    else
                        legend_label = '';
                    end
                else
                    mask = (Table.temperature == group_val) & ...
                           (Table.Total == color_val) & ...
                           (Table.rdf_type == rdf_type);
                    if ri == 1
                        legend_label = sprintf('%dH2O', color_val);
                    else
                        legend_label = '';
                    end
                end
                
                matching_idx = find(mask);
                
                if isempty(matching_idx)
                    continue
                end
                
                % Average RDF
                matching = RDF_data_filtered(matching_idx);
                [r_avg, y_avg] = average_rdf(matching, opts.ShowPMF);
                
                % Plot
                plot(r_avg, y_avg, ...
                     'Color', c, ...
                     'LineStyle', line_style, ... 
                     'LineWidth', opts. LineWidth, ...
                     'DisplayName', legend_label);
                
                if opts.Verbose
                    if opts.GroupBy == "loading"
                        fprintf('Plotting: %d H2O, %dK, %s (averaged %d)\n', ...
                                group_val, color_val, rdf_type, numel(matching));
                    else
                        fprintf('Plotting: %dK, %d H2O, %s (averaged %d)\n', ...
                                group_val, color_val, rdf_type, numel(matching));
                    end
                end
            end
        end
        
        % Set axis limits
        xlim(opts.XLim);
        
        % Apply Y-axis limits with special handling for MOF-water
        if ~isempty(opts.YLim)
            ylim(opts.YLim);
        elseif is_mof_water
            % For MOF-water interactions, limit to [0, 2]
            current_ylim = ylim;
            ylim([0, min(2, current_ylim(2))]);
        end
        
        % Add reference line
        if ~opts.ShowPMF
            yline(1, ':k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        else
            yline(0, ':k', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        
        leg = legend('Location', opts.LegendLocation);
        leg.FontSize = opts.FontSize - 1;
        
        % Save if requested
        if opts.Save
            if ~isfolder(opts. SaveDir)
                mkdir(opts.SaveDir);
            end
            
            filename = sprintf('%s.png', filename_base);
            filepath = fullfile(opts.SaveDir, filename);
            
            % Use exportgraphics for consistent high-quality output
            exportgraphics(ax, filepath, 'Resolution', 300);
            
            if opts.Verbose
                fprintf('Saved: %s\n', filepath);
            end
        end
    end
end

%% Helper: Get RDF type label for title
function label = get_rdf_type_label(rdf_types)
    % Extract the atom type from RDF type names
    % ["O-H", "O-O"] -> "O"
    % ["LinkerO-H", "LinkerO-O"] -> "LinkerO"
    % ["Mu3O-H", "Mu3O-O"] -> "Mu3O"
    
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