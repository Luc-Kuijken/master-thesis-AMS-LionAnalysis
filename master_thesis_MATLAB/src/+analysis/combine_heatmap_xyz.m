function Heatmap_data_3D = combine_heatmap_xyz(Heatmap_data, opts)
%COMBINE_HEATMAP_XYZ Combine XY, YZ, and ZX 2D slices into 3D heatmaps.
%  Heatmap_data_3D = COMBINE_HEATMAP_XYZ(Heatmap_data, Mol_values)
%  Inputs: Heatmap_data with 2D slices (XY, YZ, ZX); Mol_values table.
%  Output: Heatmap_data_3D struct with 3D interpolated fields.
%
%  Name-Value options:
%    GridSize    : Resolution of 3D grid (default: 300)
%    Method      : Interpolation method 'linear', 'natural', 'cubic' (default: 'linear')
%    Verbose     : Display progress (default: true)

    arguments
        Heatmap_data (1,:) struct
        opts.GridSize (1,1) double = 300
        opts.Method (1,1) string = "linear"
        opts.Verbose (1,1) logical = true
    end

    % Get unique combinations of H2O, H3O, T, fs
    Table = struct2table(Heatmap_data);
    XYZ_combined = unique(Table(:, {'H2O','H3O','temperature','sampling_freq'}), 'rows');
    n_combined = height(XYZ_combined);
    
    % Preallocate 3D struct
    Heatmap_data_3D(n_combined,1) = struct('name',"", 'H2O',NaN, 'H3O',NaN, ...
                                         'temperature',NaN, 'sampling_freq',NaN, ...
                                         'data3D',[]);

    for i = 1:n_combined
        H2O = XYZ_combined.H2O(i);
        H3O = XYZ_combined.H3O(i);
        T = XYZ_combined.temperature(i);
        fs = XYZ_combined.sampling_freq(i);
        
        if opts.Verbose
            fprintf("Combining 3D heatmap: %dH2O_%dH3O_%dK_%dfs\n", H2O, H3O, T, fs);
        end
        
        % Find matching 2D slices
        mask = (Table.H2O == H2O) & (Table.H3O == H3O) & ...
               (Table.temperature == T) & (Table.sampling_freq == fs);
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
            warning("Missing 2D slices for %dH2O_%dH3O_%dK_%dfs.\n" + ...
                " Skipping this heatmap.", H2O, H3O, T, fs);
            continue
        end
        
        % Extract coordinates and values
        % XY slice: Position_x, Position_y, value
        x_xy = XY.Position_x; 
        y_xy = XY.Position_y; 
        v_xy = XY.value;
        % YZ slice: Position_y, Position_z, value
        y_yz = YZ.Position_y; 
        z_yz = YZ.Position_z; 
        v_yz = YZ.value;
        % ZX slice: Position_z, Position_x, value
        z_zx = ZX.Position_z; 
        x_zx = ZX.Position_x; 
        v_zx = ZX.value;
        
        % Define regular 3D grid
        xlin = linspace(min([x_xy; x_zx]), max([x_xy; x_zx]), opts.GridSize);
        ylin = linspace(min([y_xy; y_yz]), max([y_xy; y_yz]), opts.GridSize);
        zlin = linspace(min([z_yz; z_zx]), max([z_yz; z_zx]), opts.GridSize);
        [X, Y, Z] = meshgrid(xlin, ylin, zlin);
        
        % Create scatteredInterpolant objects for each 2D slice
        F_xy = scatteredInterpolant(x_xy, y_xy, v_xy, char(opts.Method), 'none');
        F_yz = scatteredInterpolant(y_yz, z_yz, v_yz, char(opts.Method), 'none');
        F_zx = scatteredInterpolant(z_zx, x_zx, v_zx, char(opts.Method), 'none');
        
        % Interpolate each 2D slice into the 3D grid
        Vxy = NaN(size(X));
        for k = 1:numel(zlin)
            Vxy(:,:,k) = F_xy(X(:,:,k), Y(:,:,k));
        end
        Vyz = NaN(size(Y));
        for k = 1:numel(xlin)
            Vyz(:,k,:) = F_yz(squeeze(Y(:,k,:)), squeeze(Z(:,k,:)));
        end
        Vzx = NaN(size(Z));
        for k = 1:numel(ylin)
            Vzx(k,:,:) = F_zx(squeeze(Z(k,:,:)), squeeze(X(k,:,:)));
        end
        
        % Average the three interpolations
        V = mean(cat(4, Vxy, Vyz, Vzx), 4, 'omitnan');
        
        % Build 3D struct name
        struct_name = sprintf("Heatmap_%dH2O_%dH3O_%dK_%dfs_3D_XYZ", H2O, H3O, T, fs);
        
        % Fill 3D struct
        Heatmap_data_3D(i).name           = struct_name;
        Heatmap_data_3D(i).H2O            = H2O;
        Heatmap_data_3D(i).H3O            = H3O;
        Heatmap_data_3D(i).temperature    = T;
        Heatmap_data_3D(i).sampling_freq  = fs;
        Heatmap_data_3D(i).data3D         = struct('X', X, 'Y', Y, 'Z', Z, 'V', V);
    end
    
    % Remove empty entries
    valid = ~cellfun(@isempty, {Heatmap_data_3D.data3D});
    Heatmap_data_3D = Heatmap_data_3D(valid);
end