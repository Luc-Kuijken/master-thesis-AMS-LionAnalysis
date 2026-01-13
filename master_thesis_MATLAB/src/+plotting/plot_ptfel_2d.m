function fh = plot_ptfel_2d(PTFEL_data, opts)
%PLOT_PTFEL_2D Plot 2D PT-FEL heatmaps for water and MOF.
%  fh = PLOT_PTFEL_2D(PTFEL_data, Name=Value)
%
%  Creates separate 2D heatmap figures for water-water and MOF-water PT-FEL.
%  One figure per (loading, temperature) combination.
%
%  Inputs:
%    PTFEL_data : struct array from io.load_ptfel_struct
%
%  Name-Value options:
%    XLim          : x-axis limits (delta) (default: [-1.5, 1.5])
%    YLim          : y-axis limits (O-O distance) (default: [2.2, 3.5])
%    Colormap      : colormap name (default: "turbo")
%    Save          : Save figures (default: false)
%    SaveDir       : Save directory (default: "figures/ptfel")
%    Verbose       : Print messages (default: true)
%    
%    FigSize       : Figure size [width, height] (default: [700, 500])
%    FontSize      : Font size (default: 12)
%    TitleFontSize : Title font size (default: 14)

    arguments
        PTFEL_data (1,:) struct
        
        opts.XLim (1,2) double = [-1.5, 1.5]
        opts.YLim (1,2) double = [2.2, 3.5]
        opts.Colormap (1,1) string = "turbo"
        opts.Save (1,1) logical = false
        opts.SaveDir (1,1) string = "figures/ptfel"
        opts.Verbose (1,1) logical = true
        
        opts.FigSize (1,2) double = [700, 500]
        opts.FontSize (1,1) double = 12
        opts.TitleFontSize (1,1) double = 14
    end
    
    % Filter for 2D data only
    Table = struct2table(PTFEL_data);
    mask_2d = [PTFEL_data.is_2d]';
    PTFEL_2d = PTFEL_data(mask_2d);
    
    if isempty(PTFEL_2d)
        warning('No 2D PTFEL data found.');
        fh = [];
        return
    end
    
    % Plot water-water and MOF-water separately
    fh_water = plot_single_2d_ptfel(PTFEL_2d, "water", opts);
    fh_mof = plot_single_2d_ptfel(PTFEL_2d, "mof", opts);
    
    % Combine figure handles
    fh = [fh_water; fh_mof];
end

%% Helper: Plot 2D PT-FEL for a specific type
function fh = plot_single_2d_ptfel(PTFEL_2d, ptfel_type, opts)
    % Filter by type
    Table_2d = struct2table(PTFEL_2d);
    mask_type = [PTFEL_2d.type]' == ptfel_type;
    data_filtered = PTFEL_2d(mask_type);
    
    if isempty(data_filtered)
        if opts.Verbose
            warning('Could not extract 2D PT-FEL data for %s', ptfel_type);
        end
        fh = gobjects(0);
        return
    end
    
    % Get unique (loading, temperature) combinations
    Table_filtered = struct2table(data_filtered);
    unique_combos = unique(Table_filtered(: , {'H2O', 'H3O', 'temperature'}), 'rows');
    n_combos = height(unique_combos);
    
    % Preallocate figure handles
    fh = gobjects(n_combos, 1);
    
    % Loop through each combination
    for i = 1:n_combos
        H2O_val = unique_combos.H2O(i);
        H3O_val = unique_combos.H3O(i);
        T_val = unique_combos.temperature(i);
        loading = H2O_val + H3O_val;
        
        if opts.Verbose
            fprintf('Plotting 2D %s PTFEL:  %dH2O_%dH3O, %dK\n', ...
                   ptfel_type, H2O_val, H3O_val, T_val);
        end
        
        % Find matching data
        idx = find(([data_filtered.H2O] == H2O_val) & ...
                  ([data_filtered.H3O] == H3O_val) & ...
                  ([data_filtered.temperature] == T_val));
        
        if isempty(idx)
            continue
        end
        
        data = data_filtered(idx);
        
        % Extract 2D data
        delta = data.data.delta;
        distance = data.data.distance;
        fe = data.data.free_energy;
        
        % Shift free energy to zero minimum
        fe = fe - min(fe);
        
        % Create pivot table for heatmap
        % Group by (delta, distance) and average if duplicates exist
        [unique_pairs, ~, ic] = unique([delta, distance], 'rows');
        fe_grouped = accumarray(ic, fe, [], @mean);
        
        % Reshape to grid
        unique_delta = unique(delta);
        unique_dist = unique(distance);
        
        n_delta = numel(unique_delta);
        n_dist = numel(unique_dist);
        
        FE_grid = NaN(n_dist, n_delta);
        
        for k = 1:numel(fe_grouped)
            d_val = unique_pairs(k, 1);
            dist_val = unique_pairs(k, 2);
            
            [~, d_idx] = min(abs(unique_delta - d_val));
            [~, dist_idx] = min(abs(unique_dist - dist_val));
            
            FE_grid(dist_idx, d_idx) = fe_grouped(k);
        end
        
        % Create figure
        type_str = sprintf('%s- Water-Water', upper(ptfel_type(1)) + ptfel_type(2:end));
        % FIX: Ensure MOF specific naming
        if lower(ptfel_type) == "mof"
            type_str = "MOF-Water";
        end
        
        fig_name = sprintf('2D PTFEL %s - %d H2O - %d K', type_str, loading, T_val);
        fh(i) = figure('Name', fig_name, ...
                      'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]);
        ax = gca; hold on; grid on; box on;
        
        % --- PLOT HEATMAP WITH TRANSPARENCY ---
        h = imagesc(unique_delta, unique_dist, FE_grid);

        % Make both NaN AND Inf transparent (removes the "red/blue" background)
        set(h, 'AlphaData', isfinite(FE_grid)); 

        set(ax, 'Color', 'none'); % Ensure axes background is clear
        set(ax, 'YDir', 'normal');
        
        % Colormap and colorbar
        colormap(ax, opts.Colormap);
        cb = colorbar;
        cb.Label.String = 'Free Energy (kT)';
        cb.Label.FontSize = opts.FontSize;
        
        % Labels and formatting
        xlabel('\delta (Å)', 'FontSize', opts.FontSize, 'FontWeight', 'bold');
        ylabel('O \cdots O Distance (Å)', 'FontSize', opts.FontSize, 'FontWeight', 'bold');
        title(fig_name, 'FontSize', opts.TitleFontSize, 'FontWeight', 'bold');
        
        xlim(opts.XLim);
        ylim(opts.YLim);
        
        ax.FontSize = opts.FontSize;
        
        % Save if requested
        if opts.Save
            if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
            filename = sprintf('ptfel_2d_%s_%dH2O_%dH3O_%dK.png', ...
                             lower(ptfel_type), H2O_val, H3O_val, T_val);
            filepath = fullfile(opts.SaveDir, filename);
            exportgraphics(ax, filepath, 'Resolution', 300);
            if opts.Verbose, fprintf('Saved: %s\n', filepath); end
        end
    end
end