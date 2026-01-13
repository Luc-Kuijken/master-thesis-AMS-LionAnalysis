function fh = plot_2D_heatmap(Heatmap_data, opts)
%PLOT_2D_HEATMAP Plot 2D heatmap slices as scatter plots and/or 3D cube view.
%  fh = PLOT_2D_HEATMAP(Heatmap_data, Name=Value)

    arguments
        Heatmap_data (1,:) struct
        
        % Data filtering options
        opts.H2O (:,1) double = []
        opts.Temperature (:,1) double = []
        opts.SamplingFreq (:,1) double = []
        
        % Heatmap-specific options
        opts.MaxValue (1,1) double = NaN
        opts.HeatmapAlpha (1,1) double = 0
        opts.AtomAlpha (1,1) double = 0
        opts.MarkerSize (1,1) double = 3
        
        % Plot control options
        opts.PlotType (1,1) string = "XY"
        opts.DisplayMode (1,1) string = "off"
        opts.SaveDir (1,1) string = "figures/heatmaps"
        opts.MOFFile (1,1) string = "data/UiO-66.xyz"
        
        % Visibility options
        opts.ShowMOF (1,1) logical = false
        opts.ShowColorbar (1,1) logical = true
        opts.ShowLegend (1,1) logical = true
        opts.ShowTitle (1,1) logical = true
        opts.Verbose (1,1) logical = true
        opts.Save (1,1) logical = false
        
        % Style options
        opts.FontSize (1,1) double = 12
        opts.TitleFontSize (1,1) double = 14
        opts.LineWidth (1,1) double = 1.0
        opts.AxisSize (1,:) double = [380, 380]  % [width, height] in pixels
        opts.FigSize (1,:) double = [500, 420]   % [width, height] in pixels
        opts.DPI (1,1) double = 300
    end
    
    % Load MOF structure if requested
    if isfile(opts.MOFFile)
        mof = io.load_mof_structure(opts.MOFFile);
        
        if opts.Verbose
            fprintf('Loaded MOF: %s (%d atoms, unit cell: %.2f Å)\n', ...
                    mof.name, mof.n_atoms, mof.unit_cell(1,1));
        end
        
        % Plot the MOF on its own
        if opts.ShowMOF
            figure('Name', 'MOF Unit Cell - 3D');
            hold on; grid on; box on;
            axis equal;
            xlabel('X [Å]'); ylabel('Y [Å]'); zlabel('Z [Å]');
            title('UiO-66 Unit Cell');
            view(135, 35);
            axis equal padded

            plotting.plot_mof_frame(mof, ...
                Projection="3D", ...
                AtomScale=1.5);
        end
    else
        mof = [];
        warning('MOF file not found: %s. Skipping MOF overlay.', opts.MOFFile);
    end
    
    % Get unique combinations from data
    Table = struct2table(Heatmap_data);
    unique_tags = unique(Table(:, {'H2O','temperature','sampling_freq'}), 'rows');
    
    % Apply filters
    if ~isempty(opts.H2O)
        unique_tags = unique_tags(ismember(unique_tags.H2O, opts.H2O), :);
    end
    if ~isempty(opts.Temperature)
        unique_tags = unique_tags(ismember(unique_tags.temperature, opts.Temperature), :);
    end
    if ~isempty(opts.SamplingFreq)
        unique_tags = unique_tags(ismember(unique_tags.sampling_freq, opts.SamplingFreq), :);
    end
    
    n_combined = height(unique_tags);
    
    if n_combined == 0
        error('No matching heatmap data found for the specified filters.');
    end
    
    % Preallocate figure handles
    fh = gobjects(n_combined * 4, 1);
    fig_idx = 0;
    
    % Loop through all combinations
    for i = 1:n_combined
        H2O_val = unique_tags.H2O(i);
        T_val = unique_tags.temperature(i);
        fs_val = unique_tags.sampling_freq(i);
        
        if opts.Verbose
            fprintf('Plotting: %dH2O, %dK, %dfs\n', H2O_val, T_val, fs_val); 
        end
        
        % Find matching slices for this combination
        mask = (Table.H2O == H2O_val) & ...
               (Table.temperature == T_val) & ...
               (Table.sampling_freq == fs_val);
        slices = Heatmap_data(mask);
        
        % Identify XY, YZ, ZX slices
        XY = []; YZ = []; ZX = [];
        for ii = 1:numel(slices)
            coords = sort([char(slices(ii).Position_1), char(slices(ii).Position_2)]);
            switch string(coords)
                case "XY"
                    XY = slices(ii).data;
                case "YZ"
                    YZ = slices(ii).data;
                case "XZ"
                    ZX = slices(ii).data;
            end
        end
    
        % Check if all three slices are present
        if isempty(XY) || isempty(YZ) || isempty(ZX)
            warning('Missing 2D slices for %dH2O_%dK_%dfs. Skipping.', H2O_val, T_val, fs_val);
            continue
        end
    
        % Extract data
        a(:,1) = XY.Position_x; b(:,1) = XY.Position_y; v(:,1) = XY.value;
        a(:,2) = YZ.Position_y; b(:,2) = YZ.Position_z; v(:,2) = YZ.value;
        a(:,3) = ZX.Position_z; b(:,3) = ZX.Position_x; v(:,3) = ZX.value;

        % Normalize density to [0, 1] (relative density)
        if ~isnan(opts.MaxValue)
            % Use MaxValue as the reference maximum
            v_raw_max = opts.MaxValue;
        else
            % Use the actual maximum from all slices
            v_raw_max = max([v(:,1); v(:,2); v(:,3)]);
        end
        
        % Normalize each slice
        v(:,1) = v(:,1) / v_raw_max;
        v(:,2) = v(:,2) / v_raw_max;
        v(:,3) = v(:,3) / v_raw_max;
        
        % Clamp values to [0, 1]
        v = max(0, min(1, v));
        
        % Set normalized limits
        v_min = 0;
        v_max = 1;
    
        % Plot individual slices
        slice_names = ["XY", "YZ", "ZX"];
        for iii = 1:length(slice_names)
            slice = slice_names(iii);
            if ismember(opts.PlotType, [slice, "XYZ", "all"])
                fig_idx = fig_idx + 1;
                fh(fig_idx) = plot_single_slice(a(:,iii), b(:,iii), ...
                    v(:,iii), slice, opts, v_min, v_max, mof, slices(1)); 
            end
        end
        
        % Plot 3D cube view
        if ismember(opts.PlotType, ["3D", "all"])
            fig_idx = fig_idx + 1;
            fh(fig_idx) = plot_cube_view(a(:,1), b(:,1), v(:,1), ...
                                         a(:,2), b(:,2), v(:,2), ...
                                         a(:,3), b(:,3), v(:,3), ...
                                         opts, v_min, v_max, mof, slices(1));
        end
        
        % Return error when not plotted
        if ~ismember(opts.PlotType, ["XY", "YZ", "ZX", "XYZ", "3D", "all"])
            error('Choose PlotType from: ["XY", "YZ", "ZX", "XYZ", "3D", "all"]')
        end
    end
    
    % Remove empty figure handles
    fh = fh(1:fig_idx);
end

%% Helper function: Apply consistent styling to axes
function apply_axis_style(ax, opts)
    ax.FontSize = opts.FontSize;
    ax.LineWidth = opts.LineWidth;
    ax.XLabel.FontSize = opts.TitleFontSize;
    ax.YLabel.FontSize = opts.TitleFontSize;
    if isprop(ax, 'ZLabel')
        ax.ZLabel.FontSize = opts.TitleFontSize;
    end
    ax.Title.FontSize = opts.TitleFontSize;
end

%% Helper function: Set fixed axis size (ensures consistent plot area)
function set_fixed_axis_size(fh, ax, axis_size, has_colorbar)
    if isempty(axis_size) || numel(axis_size) < 2
        return  % Don't modify if no fixed size specified
    end
    
    axis_width = axis_size(1);
    axis_height = axis_size(2);
    
    % Margins (in pixels)
    margin_left = 75;
    margin_right = 25;
    margin_bottom = 65;
    margin_top = 55;
    colorbar_width = 85;
    
    % Calculate figure size
    if has_colorbar
        fig_width = margin_left + axis_width + colorbar_width + margin_right;
    else
        fig_width = margin_left + axis_width + margin_right;
    end
    fig_height = margin_bottom + axis_height + margin_top;
    
    % Set figure size
    fh.Position(3:4) = [fig_width, fig_height];
    
    % Set axis position in pixels
    ax.Units = 'pixels';
    ax.Position = [margin_left, margin_bottom, axis_width, axis_height];
end

%% Helper function: plot single 2D slice
function fh = plot_single_slice(x, y, v, slice_name, opts, v_min, v_max, mof, metadata)

    % Generate filename and figname
    [filename, fig_name] = generate_filename(metadata, slice_name, ("_2D_" + slice_name));
    filepath = fullfile(opts.SaveDir, filename);
    
    % Check if image mode and file exists
    if opts.DisplayMode == "image" && isfile(filepath)
        fh = display_png(filepath, fig_name);
        return
    elseif opts.DisplayMode == "interactive"
        visibility = 'on';
    else
        visibility = 'off';
    end
    
    % Create figure
    fh = figure('Name', fig_name, 'Visible', visibility, ... 
                'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]);
    ax = gca; hold on; grid on; box on;
    
    % Scatter plot with color and optionally alpha
    s = scatter(x, y, opts.MarkerSize, v, 'filled');
    if opts.HeatmapAlpha < 1
        alpha_vals = calc_alpha(v, v_min, v_max, opts.HeatmapAlpha);
        s.AlphaData = alpha_vals;
        s.MarkerFaceAlpha = 'flat';
    end
    
    % Colormap and limits
    colormap(ax, "turbo");
    caxis([v_min, v_max]);
    
    % Add colorbar if requested
    if opts.ShowColorbar
        cb = colorbar;
        cb.Label.String = 'Relative Density';
        cb.Label.FontSize = opts.TitleFontSize;
        cb.FontSize = opts.FontSize;
    end
    
    % Set fixed axis size AFTER adding colorbar
    set_fixed_axis_size(ax, opts, opts.ShowColorbar);
    
    % Labels
    switch slice_name
        case 'XY'
            xlabel('Position X [Å]');
            ylabel('Position Y [Å]');
        case 'YZ'
            xlabel('Position Y [Å]');
            ylabel('Position Z [Å]');
        case 'ZX'
            xlabel('Position Z [Å]');
            ylabel('Position X [Å]');
    end
    
    if opts.ShowTitle
        title(fig_name);
    end
    axis equal padded;
    
    % Apply consistent styling
    apply_axis_style(ax, opts);
    
    % Overlay MOF unit cell with or without MOF atoms
    if ~isempty(mof)
        plotting.plot_mof_frame(mof, ...
            Projection=slice_name, ...
            LineStyle="--", ...
            AtomScale=1.5, ...
            AtomAlpha=opts.AtomAlpha, ...
            LineWidth=opts.LineWidth);
    end
    
    % Save if requested
    if opts.Save || opts.DisplayMode == "image"
        if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
        exportgraphics(fh, filepath, 'Resolution', opts.DPI);
        if opts.DisplayMode ~= "interactive", close(fh); end
        if opts.Verbose, fprintf('Saved: %s\n', filepath); end
        
        % If image mode, close and reopen with image
        if opts.DisplayMode == "image"
            fh = display_png(filepath, fig_name);
        end
    end
end

%% Helper function: plot 3D cube view with three faces
function fh = plot_cube_view(x_xy, y_xy, v_xy, y_yz, z_yz, v_yz, ...
                             z_zx, x_zx, v_zx, opts, v_min, v_max, mof, metadata)
    
    % Generate filename using helper
    [filename, fig_name] = generate_filename(metadata, '', "_3D_XYZ");
    filepath = fullfile(opts.SaveDir, filename);
    
    % Check if image mode and file exists
    if opts.DisplayMode == "image" && isfile(filepath)
        fh = display_png(filepath, fig_name);
        return
    elseif opts.DisplayMode == "interactive"
        visibility = 'on';
    else
        visibility = 'off';
    end
    
    % Create figure with specified size
    fh = figure('Name', fig_name, 'Visible', visibility, ... 
                'Position', [100, 100, opts.FigSize(1), opts.FigSize(2)]);
    ax = gca; hold on; grid on; box on;
    
    % Get data bounds
    x_max = max([x_xy; x_zx]);
    y_max = max([y_xy; y_yz]);
    z_max = max([z_yz; z_zx]);
    
    % XY slice (top face at max Z)
    s1 = scatter3(x_xy, y_xy, repmat(z_max, size(x_xy)), opts.MarkerSize, v_xy, 'filled');
    % YZ slice (right face at max X)
    s2 = scatter3(repmat(x_max, size(y_yz)), y_yz, z_yz, opts.MarkerSize, v_yz, 'filled');
    % ZX slice (back face at max Y)
    s3 = scatter3(x_zx, repmat(y_max, size(z_zx)), z_zx, opts.MarkerSize, v_zx, 'filled');
    
    % Optional: add alpha
    if opts.HeatmapAlpha < 1
        factor_3d = 0.4; % this factor is to correct the opacity for the diagonal view
        alpha_xy = calc_alpha(v_xy, v_min, v_max, opts.HeatmapAlpha * factor_3d);
        s1.AlphaData = alpha_xy;    s1.MarkerFaceAlpha = 'flat';
        alpha_yz = calc_alpha(v_yz, v_min, v_max, opts.HeatmapAlpha * factor_3d);
        s2.AlphaData = alpha_yz;    s2.MarkerFaceAlpha = 'flat';
        alpha_zx = calc_alpha(v_zx, v_min, v_max, opts.HeatmapAlpha * factor_3d);
        s3.AlphaData = alpha_zx;    s3.MarkerFaceAlpha = 'flat';
    end
    
    % Colormap and limits
    colormap(ax, "turbo");
    caxis([v_min, v_max]);
    
    if opts.ShowColorbar
        cb = colorbar;
        cb.Label.String = 'Relative Density';
        cb.Label.FontSize = opts.TitleFontSize;
        cb.FontSize = opts.FontSize;
    end
    
    % Labels
    xlabel('Position X [Å]');
    ylabel('Position Y [Å]');
    zlabel('Position Z [Å]');
    ax.XLabel.Rotation = 30;
    ax.YLabel. Rotation = -30;
    
    if opts.ShowTitle
        title(fig_name);
    end
    
    % Add padding and set view
    axis equal tight;
    pad = 0.07;
    xlim([0, x_max*(1+pad)]);
    ylim([0, y_max*(1+pad)]);
    zlim([0, z_max*(1+pad)]);
    view(135, 35);  
    
    % Apply consistent styling
    apply_axis_style(ax, opts);
    
    % Overlay 3D MOF unit cell without MOF atoms for 3D heatmap
    if ~isempty(mof)
        plotting.plot_mof_frame(mof, ...
            Projection="3D", ...
            LineStyle="--", ...
            AtomAlpha=0, ...
            LineWidth=opts.LineWidth);
    end
    
    % Save if requested
    if opts.Save || opts.DisplayMode == "image"
        if ~isfolder(opts.SaveDir), mkdir(opts.SaveDir); end
        exportgraphics(fh, filepath, 'Resolution', opts.DPI);
        if opts.DisplayMode ~= "interactive", close(fh); end
        if opts.Verbose, fprintf('Saved: %s\n', filepath); end
        
        % If image mode, close and reopen with image
        if opts.DisplayMode == "image"
            fh = display_png(filepath, fig_name);
        end
    end
end

%% Extra Helper functions: 

% generate filename and get time with corresponding unit
function [filename, fig_name] = generate_filename(metadata, slice_name, file_suffix)
    % Calculate time and unit
    time = metadata.sampling_freq * 500;
    if time < 1e5
        time = time/1e3; unit = 'ps';
    else
        time = time/1e6; unit = 'ns';
    end
    
    % Generate filename
    filename = sprintf('AvOPs_%dH2O_0H3O_%dK_%d%s%s.png', ...
                      metadata.H2O, metadata.temperature, time, unit, file_suffix);
    
    % Generate figure name
    if nargin > 1 && ~isempty(slice_name)
        fig_name = sprintf('2D %s slice - %dH2O 0H3O %dK %d%s', ...
                          slice_name, metadata.H2O, metadata.temperature, time, unit);
    else
        fig_name = sprintf('3D Heatmap - %dH2O 0H3O %dK %d%s', ...
                          metadata.H2O, metadata.temperature, time, unit);
    end
end

% calculate alpha values
function alpha = calc_alpha(v, v_min, v_max, alpha_value)
    % Normalize to [0, 1] and apply power transform
    alpha = (v - v_min) / (v_max - v_min);
    alpha = max(0, min(1, alpha));
    power = 1-alpha_value; % invert because alpha = 1 is achieved by power = 0
    alpha = alpha.^power;
end

% display PNG image
function fh = display_png(filepath, fig_name)
    % Fast display of pre-rendered PNG
    fh = figure('Name', fig_name);
    imshow(imread(filepath));
    title(fig_name);
end