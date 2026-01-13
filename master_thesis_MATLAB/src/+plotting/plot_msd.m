function fh = plot_msd(Diffusion_data, Mol_values, opts)
%PLOT_MSD Plot MSD vs time (and linear fits) per H2O/H3O group.
%  fh = PLOT_MSD(Diffusion_data, Mol_values, Name=Value)
%
%  This function groups simulations by the (H2O,H3O) combinations listed
%  in Mol_values and produces one figure per combination showing:
%    - MSD vs time for each simulation (different temperatures)
%    - The linear fit used to estimate diffusion (if available)
%
%  Inputs:
%    Diffusion_data : 1byN struct array with fields:
%                      .name, .H2O, .H3O, .temperature, .data.time,
%                      .data.msd_a, .fit_coeffs (slope/intercept)
%    Mol_values     : table or Nx2 table containing columns H2O and H3O
%
%  Name-Value options (opts):
%    Scale          : "linear" (default) or "log" - applied to both axes
%    Save           : logical false (default) - save generated figures
%    SaveDir        : directory to write saved figures ("figures")
%    LegendLocation : location passed to legend (default "northwest")
%    FigSize        : Figure size [width, height] (default:  [800, 400])
%    FontSize       : Font size for labels (default: 14)
%    TitleFontSize  : Font size for title (default: 16)
%    LineWidth      : Line width (default: 1.5)
%
%  Output:
%    fh : array of figure handles

    arguments
        Diffusion_data (1,:) struct
        Mol_values (:,3) table
        opts.Scale (1,1) string = "linear"
        opts.Save (1,1) logical = false
        opts.SaveDir (1,1) string = "figures"
        opts.LegendLocation (1,1) string = "northwest"
        opts.FigSize (1,2) double = [800, 400]
        opts.FontSize (1,1) double = 16
        opts.TitleFontSize (1,1) double = 18
        opts.LineWidth (1,1) double = 1.5
    end

    % Number of unique H2O/H3O combinations to plot
    num_sims = numel(Mol_values.H2O);
    fh = gobjects(num_sims,1);

    for i = 1:num_sims
        % Select simulations that match the current H2O value
        idx   = [Diffusion_data.H2O] == Mol_values.H2O(i);
        group = Diffusion_data(idx);

        % Create figure
        fh(i) = figure('Name', sprintf('MSD %s scale, %gH2O %gH3O', ...
            opts.Scale, Mol_values{i,1:2}), ...
            'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]);

        % Prepare axes and plotting defaults
        ax = gca; hold on; grid on; box on
        ColOrd = ax.ColorOrder;
        % Set Scale and axis
        if opts.Scale == "log"
            set(ax,'XScale', opts.Scale, 'YScale', opts.Scale)
            axis([1e-1 1e3 1e-1 1e3])
        else
            if opts.Scale ~= "linear"
            	warning('Invalid Scale value. Using "linear".');
                opts.Scale = 'linear';
            end
            set(ax,'XScale', opts.Scale, 'YScale', opts.Scale)
            axis([0 1000 0 600])
        end
        xlabel('time [ps]', 'FontSize', opts.FontSize); 
        ylabel('MSD [Å^2]', 'FontSize', opts.FontSize);
        title(sprintf('MSD %s scale, %gH2O %gH3O', ... 
            opts.Scale, Mol_values{i,1:2}), ...
            'FontSize', opts.TitleFontSize, 'FontWeight', 'bold');
        ax.FontSize = opts.FontSize;

        % Loop over all simulations in this group and plot MSD + fit
        for ii = 1:numel(group)
            % Cycle through color order
            c = ColOrd(mod(ii-1, size(ColOrd,1)) + 1, :);

            % Extract time and MSD arrays for plotting
            t = group(ii).data.time;
            y = group(ii).data.msd_a;

            % Plot MSD data
            plot(t, y, 'Color', c, 'LineWidth', opts.LineWidth, ...
                'DisplayName', sprintf('%g K', group(ii).temperature));

            % Plot the linear fit
            if isfield(group(ii), 'fit_coeffs') && ~isempty(group(ii).fit_coeffs)
                % fit_coeffs is expected as [slope, intercept] from polyfit
                yfit = t .* group(ii).fit_coeffs(1) + group(ii).fit_coeffs(2);
                plot(t, yfit, '--', 'Color', c, 'HandleVisibility','off', ... 
                     'LineWidth', opts. LineWidth * 0.8);
            end
        end

        % Show legend
        leg = legend('Location', opts.LegendLocation);
        leg.FontSize = opts.FontSize - 1;
        % Save as a PNG image using exportgraphics
        if opts.Save
            if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
            out = fullfile(opts.SaveDir, sprintf('msd_%s_%g_H2O_%g_H3O.png', ...
                opts.Scale, Mol_values{i,1:2}));
            exportgraphics(ax, out)
        end
    end
end


