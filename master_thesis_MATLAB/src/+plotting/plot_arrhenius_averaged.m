function fh = plot_arrhenius_averaged(Arrhenius_params, k_B, opts)
%PLOT_ARRHENIUS_AVERAGED Plot averaged Arrhenius data with RMSE uncertainty bands.
%  fh = PLOT_ARRHENIUS_AVERAGED(Arrhenius_params, k_B, Name=Value)
%  Uses averaged Arrhenius parameters and RMSE for uncertainty quantification.
%  Name-Value:
%    ShowFitError  : Show RMSE uncertainty bands around fit (default:  true)
%    ErrorAlpha    : Transparency for RMSE bands (default: 0.2)
%    MarkerSize    : Marker size (default: 6)
%    LineWidth     : Line width (default: 1.5)
%    FigSize       : Figure size [width, height] (default: [800, 400])
%    FontSize      : Font size for labels (default:  14)
%    TitleFontSize : Font size for title (default: 16)

    arguments
        Arrhenius_params (:,7) table
        k_B (1,1) double
        opts.Save (1,1) logical = false
        opts.SaveDir (1,1) string = "figures"
        opts.LegendLocation (1,1) string = "northeast"
        opts.MarkerSize (1,1) double = 6
        opts.LineWidth (1,1) double = 1.5
        opts.ShowFitError (1,1) logical = true
        opts.ErrorAlpha (1,1) double = 0.2
        opts.FigSize (1,2) double = [800, 400]
        opts.FontSize (1,1) double = 16
        opts.TitleFontSize (1,1) double = 18
    end
    
    % Create figure
    fh = figure('Name', 'Arrhenius plot - Averaged per loading', ...
    	'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]);
    ax = gca; hold on; grid on; box on;
    xlabel('1/T [K^{-1}]', 'FontSize', opts.FontSize); 
    ylabel('ln(D) [ln(m^2 s^{-1})]', 'FontSize', opts. FontSize);
    title('Arrhenius fit: Averaged per loading', ...
          'FontSize', opts.TitleFontSize, 'FontWeight', 'bold');
    ax.FontSize = opts.FontSize;
    
    ColOrd = [ax.ColorOrder(4:end,:); [1.0, 0.6, 0.0]];
    
    % Get number of loadings
    n_loadings = height(Arrhenius_params);
    
    % Plot each loading
    for i = 1:n_loadings
        loading_val = Arrhenius_params.Total(i);
        
        % Get stored averaged points and RMSE
        T_avg = Arrhenius_params.T_avg{i};
        D_avg = Arrhenius_params.D_avg{i};
        RMSE = Arrhenius_params.RMSE(i);
        
        % Convert to Arrhenius coordinates
        X = 1 ./ T_avg;
        Y = log(D_avg);
        
        % Select color
        c = ColOrd(mod(i-1, size(ColOrd,1)) + 1, :);
        
        % Plot fit line using parameters
        Ea = Arrhenius_params.Ea(i);
        D0 = Arrhenius_params.D0(i);
        
        Xfit = linspace(min(X)*0.95, max(X)*1.05, 100);
        Yfit = log(D0) - (Ea / k_B) * Xfit;
        
        % Upper and lower bounds with RMSE
        Yfit_upper = Yfit + RMSE;
        Yfit_lower = Yfit - RMSE;
        
        % Fill RMSE uncertainty band
        fill([Xfit, fliplr(Xfit)], [Yfit_upper, fliplr(Yfit_lower)], ...
             c, 'FaceAlpha', opts.ErrorAlpha, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
         
        % Plot fit line (with legend entry)
        plot(Xfit, Yfit, '--', 'Color', c, 'LineWidth', opts.LineWidth, ...
             'DisplayName', sprintf('%g H2O/Cell', loading_val));
        
        % Plot averaged data points (no legend entry to avoid duplicates)
        plot(X, Y, 'x', 'MarkerSize', opts.MarkerSize, 'LineWidth', ...
            opts.LineWidth, 'Color', c, 'HandleVisibility', 'off');
    end
    
    leg = legend('Location', opts.LegendLocation);
    leg.FontSize = opts.FontSize - 2;
    xlim([min(Xfit) max(Xfit)])
    
    % Save if requested
    if opts.Save
        if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
        out = fullfile(opts.SaveDir, 'arrhenius_averaged.png');
        exportgraphics(ax, out, 'Resolution', 300);
    end
end