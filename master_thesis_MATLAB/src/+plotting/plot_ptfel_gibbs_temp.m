function fh = plot_ptfel_gibbs_temp(PTFEL_params, opts)
%PLOT_PTFEL_GIBBS_TEMP Plot ΔG‡ vs temperature with linear fits per loading.
%  fh = PLOT_PTFEL_GIBBS_TEMP(PTFEL_params, Name=Value)
%
%  Shows Gibbs free energy barriers as a function of temperature.
%  Each loading gets a unique color (matching RDF color scheme).
%  Points are marked with circles, fits with lines.
%
%  Inputs: 
%    PTFEL_params : table with columns:
%                   loading, type, direction, dG_300K_kJmol, dG_350K_kJmol, 
%                   dG_400K_kJmol, dH_kJmol, dS_JmolK
%
%  Name-Value options:
%    Type          : "water" | "mof" | "both" (default: "water")
%    Save          : Save figure (default: false)
%    SaveDir       : Save directory (default: "figures/ptfel")
%    Verbose       : Print messages (default: true)
%    
%    FigSize       :  Figure size [width, height] (default: [800, 400])
%    FontSize      :  Font size (default: 14)
%    TitleFontSize : Title font size (default: 16)
%    LineWidth     : Line width (default: 1.5)
%    MarkerSize    : Marker size (default: 8)

    arguments
        PTFEL_params (:,13) table
        
        opts.Type (1,1) string = "water"
        opts.Save (1,1) logical = false
        opts.SaveDir (1,1) string = "figures/ptfel"
        opts.Verbose (1,1) logical = true
        
        opts.FigSize (1,2) double = [800, 400]
        opts.FontSize (1,1) double = 14
        opts.TitleFontSize (1,1) double = 16
        opts.LineWidth (1,1) double = 1.5
        opts.MarkerSize (1,1) double = 8
    end
    
    % Filter by type
    if opts.Type ~= "both"
        mask = PTFEL_params.type == opts.Type;
        PTFEL_params = PTFEL_params(mask, :);
    end
    
    if isempty(PTFEL_params)
        warning('No PTFEL parameters found for type=%s', opts.Type);
        fh = [];
        return
    end
    
    % Get unique loadings
    unique_loadings = unique(PTFEL_params.loading);
    
    % Use RDF color scheme (same as plot_rdf_combined)
    ax_temp = axes('Visible', 'off');
    default_colors = ax_temp.ColorOrder;
    close(gcf);
    loading_colors = [default_colors(4: end,: ); [1.0, 0.6, 0.0]];
    
    % Temperature points
    temps = [300, 350, 400];
    
    % Create figure
    fig_name = 'Gibbs Free Energy of Proton Transfer vs Temperature';
    fh = figure('Name', fig_name, ...
               'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]);
    ax = gca; hold on; grid on; box on;
    
    % Plot each loading
    for i = 1:numel(unique_loadings)
        loading = unique_loadings(i);
        
        % Select color
        c = loading_colors(mod(i-1, size(loading_colors, 1)) + 1, :);
        
        % Get data for this loading
        mask_loading = PTFEL_params.loading == loading;
        params_loading = PTFEL_params(mask_loading, :);
        
        % Average over directions if MOF (left/right)
        dG_300K = mean(params_loading.dG_300K_kJmol);
        dG_350K = mean(params_loading.dG_350K_kJmol);
        dG_400K = mean(params_loading.dG_400K_kJmol);
        
        dG_vals = [dG_300K, dG_350K, dG_400K];
        
        % Plot data points with markers
        plot(temps, dG_vals, 'o', 'Color', c, ...
             'MarkerSize', opts.MarkerSize, 'MarkerFaceColor', c, ...
             'LineWidth', 1.5, 'HandleVisibility', 'off');
        
        % Fit line:  ΔG = ΔH - T·ΔS
        dH = mean(params_loading.dH_kJmol);
        dS = mean(params_loading.dS_JmolK);
        
        % Generate fit line
        T_fit = linspace(280, 420, 100);
        dG_fit = dH - T_fit .* dS / 1000;  % Convert J/mol/K to kJ/mol
        
        plot(T_fit, dG_fit, '-', 'Color', c, 'LineWidth', opts.LineWidth, ...
             'DisplayName', sprintf('%d H_2O/Cell (\\Delta G^{\\ddagger}=%.1f kT)', loading, dG_350K));
    end
    
    % Labels and formatting
    xlabel('Temperature (K)', 'FontSize', opts.FontSize, 'FontWeight', 'bold');
    ylabel('\Delta G^{\ddagger} (kJ/mol)', 'FontSize', opts.FontSize, 'FontWeight', 'bold');
    title(fig_name, 'FontSize', opts.TitleFontSize, 'FontWeight', 'bold');
    
    xlim([280, 420]);
    
    ax.FontSize = opts.FontSize;
    legend('Location', 'northwest', 'FontSize', opts.FontSize - 2);
    
    % Save if requested
    if opts.Save
        if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
        filename = sprintf('ptfel_gibbs_vs_temp_%s.png', opts.Type);
        filepath = fullfile(opts.SaveDir, filename);
        exportgraphics(ax, filepath, 'Resolution', 300);
        if opts.Verbose, fprintf('Saved: %s\n', filepath); end
    end
end