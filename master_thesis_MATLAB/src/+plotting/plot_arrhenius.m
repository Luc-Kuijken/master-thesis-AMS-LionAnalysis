function fh = plot_arrhenius(Diffusion_data, Mol_values, k_B, opts)
%PLOT_ARRHENIUS Plot ln(D) vs 1/T with fitted lines per H2O group.
%  fh = PLOT_ARRHENIUS(Diffusion_data, H2O_values, k_B, Name=Value)
%  Name-Value:  Save (false), SaveDir ("figures"), LegendLocation ("northeast"),
%              MarkerSize (6), LineWidth (1.5), FigSize ([800,400]),
%              FontSize (14), TitleFontSize (16)

    arguments
        Diffusion_data (1,:) struct
        Mol_values (:,3) table
        k_B (1,1) double
        opts.Save (1,1) logical = false
        opts.SaveDir (1,1) string = "figures"
        opts.LegendLocation (1,1) string = "northeast"
        opts.MarkerSize (1,1) double = 6
        opts.LineWidth (1,1) double = 1.5
        opts.FigSize (1,2) double = [800, 400]
        opts.FontSize (1,1) double = 16
        opts.TitleFontSize (1,1) double = 18
    end
    
    unique_loadings = unique(Mol_values.Total);
    num_loadings = numel(unique_loadings);
    
    fh = gobjects(num_loadings, 1);

    for i = 1:num_loadings
        current_loading = unique_loadings(i);
        
        % Create new figure for this loading
        fig_name = sprintf('Arrhenius fit - loading: %g molecules', current_loading);
        fh(i) = figure('Name', fig_name, ...
        	'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]);
        ax = gca; hold on; grid on; box on;
        xlabel('1/T [K^{-1}]', 'FontSize', opts.FontSize); 
        ylabel('ln(D) [ln(m^2 s^{-1})]', 'FontSize', opts.FontSize);
        title(fig_name, 'FontSize', opts.TitleFontSize, 'FontWeight', 'bold');
        ax.FontSize = opts.FontSize;
        
        ColOrd = [ax.ColorOrder(4:end,:); [1.0, 0.6, 0.0]];
        
        % Find all compositions with this total loading
        loading_compositions = Mol_values(Mol_values.Total == current_loading, :);
        
        % Plot each composition with this loading
        for j = 1:height(loading_compositions)
            H2O_val = loading_compositions.H2O(j);
            H3O_val = loading_compositions.H3O(j);
            
            % Find all data matching this H2O and H3O combination
            idx = ([Diffusion_data.H2O] == H2O_val) & ...
                  ([Diffusion_data.H3O] == H3O_val);
            
            if ~any(idx), continue, end
            
            T = [Diffusion_data(idx).temperature];
            D = [Diffusion_data(idx).D];
            
            X = 1 ./ T;
            Y = log(D);
            
            % Select color from color order
            c = ColOrd(mod(j-1, size(ColOrd,1)) + 1, :);
            
            % Plot data points
            plot(X, Y, 'x', 'Color', c, 'MarkerSize', opts.MarkerSize, ...
                 'LineWidth', opts.LineWidth, 'DisplayName', ...
                 sprintf('%gH2O %gH3O', H2O_val, H3O_val));
            
            % Plot fit line using parameters
            kAny = find(idx, 1, 'first');
            Ea = Diffusion_data(kAny).Ea;
            D0 = Diffusion_data(kAny).D0;
            
            Xfit = linspace(min(X)*0.95, max(X)*1.05, 100);
            Yfit = log(D0) - (Ea / k_B) * Xfit;
            plot(Xfit, Yfit, '--', 'Color', c, 'LineWidth', ...
                opts.LineWidth * 0.8, 'HandleVisibility', 'off');
        end
        
        leg = legend('Location', opts.LegendLocation);
        leg.FontSize = opts.FontSize - 1;
        xlim([min(Xfit) max(Xfit)])
        
        if opts.Save
            if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
            out = fullfile(opts.SaveDir, sprintf('arrhenius_loading_%g.png', ...
            	current_loading));
            exportgraphics(ax, out);
        end
    end
end


