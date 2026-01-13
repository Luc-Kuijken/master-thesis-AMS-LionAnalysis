function fh = plot_beta(Diffusion_data, Mol_values, opts)
%PLOT_BETA Plot beta (MSD exponent) vs time grouped by loading.
%  fh = PLOT_BETA(Diffusion_data, Mol_values, Name=Value)
%
%  Beta is the time-derivative exponent of MSD, calculated as:
%    beta = d(ln(MSD))/d(ln(t))
%  For diffusive behavior, beta should approach 1.
%
%  Inputs:
%    Diffusion_data : 1byN struct array with fields:
%                      .H2O, .H3O, .temperature, .data.time, .data.beta
%    Mol_values     : table with columns H2O, H3O, and Total (loading)
%
%  Name-Value pairs (opts):
%    Scale          : "linear" (default) or "log" - applied to x-axis (time)
%    Save           : logical false (default) - save generated figures
%    SaveDir        : directory to write saved figures ("figures")
%    LegendLocation : legend placement (default "northwest")
%    LineAlpha      : transparency for lines (default 0.4)
%    LineWidth      : line width (default 0.5)
%    AvgLineWidth   : average line width (default 1.5)
%    FigSize        : Figure size [width, height] (default: [800, 400])
%    FontSize       : Font size for labels (default: 16)
%    TitleFontSize  : Font size for title (default: 18)
%
%  Output:
%    fh : array of figure handles (one per unique loading)

    arguments
        Diffusion_data (1,:) struct
        Mol_values (:,3) table
        opts.Scale (1,1) string = "linear"
        opts.Save (1,1) logical = false
        opts.SaveDir (1,1) string = "figures"
        opts.LegendLocation (1,1) string = "northwest"
        opts.LineAlpha (1,1) double = 0.4
        opts.LineWidth (1,1) double = 0.5
        opts.AvgLineWidth (1,1) double = 1.5
        opts.FigSize (1,2) double = [800, 400]
        opts.FontSize (1,1) double = 16
        opts.TitleFontSize (1,1) double = 18
    end
    
    % Get unique loadings and temperatures
    unique_loadings = unique(Mol_values.Total);
    all_temps = unique([Diffusion_data.temperature]);
    n_temps = numel(all_temps);
    
    % Preallocate figure handles
    n_loadings = numel(unique_loadings);
    fh = gobjects(n_loadings, 1);
    
    % Get the max row counts
    max_len = max(arrayfun(@(x) height(x.data), Diffusion_data));

    % Find index of the first entry with max_len rows
    idx_max = find(arrayfun(@(x) height(x.data), Diffusion_data) == max_len, 1, 'first');

    % Get the time array from that entry
    beta_time = Diffusion_data(idx_max).data.time;
    
    % Loop over each loading
    for i = 1:n_loadings
        loading_val = unique_loadings(i);
        
        % Find all H2O/H3O combinations for this loading
        loading_mask = Mol_values.Total == loading_val;
        loading_combos = Mol_values(loading_mask, :);
        
        % Create figure for this loading
        fig_name = sprintf('Beta - %g H2O', loading_val);
        fh(i) = figure('Name', fig_name, ...
            'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]);
        ax = gca; hold on; grid on; box on;
        set(ax, 'XScale', opts.Scale);
        xlabel('time [ps]', 'FontSize', opts.FontSize); 
        ylabel('\beta', 'FontSize', opts.FontSize+2, 'FontWeight', 'bold');
        title(fig_name, 'FontSize', opts.TitleFontSize, 'FontWeight', 'bold');
        ax. FontSize = opts.FontSize;
        
        % Get default color order
        ColOrd = ax.ColorOrder;
        
        % Add reference line at beta = 1 (diffusive regime)
        yline(1, '--k', 'LineWidth', 1.5, 'DisplayName', 'Diffusive (\beta=1)', ...
              'HandleVisibility', 'on');
        
        % Loop over all temperatures H2O/H3O combinations for this loading
        for ti = 1:n_temps
            temp_val = all_temps(ti);
            c = ColOrd(mod(ti-1, size(ColOrd,1)) + 1, :);
            
            % Initialize beta for this temperature
            beta_data = NaN(max_len, height(loading_combos));
            
            % Loop over all H2O/H3O combinations for this loading
            for j = 1:height(loading_combos)
                H2O_val = loading_combos.H2O(j);
                H3O_val = loading_combos.H3O(j);

                % Find simulation matching this temp and H2O/H3O combo
                idx = ([Diffusion_data.H2O] == H2O_val) & ...
                      ([Diffusion_data.H3O] == H3O_val) & ...
                      ([Diffusion_data.temperature] == temp_val);
            
                % Plot and accumulate for each simulation at this temperature
                t = Diffusion_data(idx).data.time;
                beta = Diffusion_data(idx).data.beta;

                % Plot individual trace with transparency
                p = plot(t, beta, ':', 'Color', c, 'LineWidth', ...
                    opts.LineWidth, 'HandleVisibility', 'off');
                p.Color(4) = opts.LineAlpha;
                
                beta_len = length(beta);
                beta_data(1:beta_len,j) = beta(:);
            end
            
            beta_avg = mean(beta_data, 2, 'omitnan');
            plot(beta_time, beta_avg, '-', 'Color', c, 'LineWidth', opts.AvgLineWidth, ...
                 'DisplayName', sprintf('%g K (avg)', temp_val));
        end
        
        % Show legend and set axis
        leg = legend('Location', opts.LegendLocation);
        leg.FontSize = opts.FontSize - 1;
        axis([0 1e3 0 2]);
        
        % Save if requested
        if opts.Save
            if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
            out = fullfile(opts.SaveDir, sprintf('beta_loading_%g.png', loading_val));
            exportgraphics(ax, out, 'Resolution', 300);
        end
    end
end