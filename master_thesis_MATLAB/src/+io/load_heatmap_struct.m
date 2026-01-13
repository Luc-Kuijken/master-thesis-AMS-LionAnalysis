function [Heatmap_data, Mol_values] = load_heatmap_struct(folderPath, opts)
%LOAD_HEATMAP_STRUCT Scan a folder of heatmap .dat files and build Heatmap_data.
%  [Heatmap_data, Mol_values] = io.load_heatmap_struct("data/raw", Pattern="*_2D_*.dat", Verbose=true)

    arguments
        folderPath (1,1) string
        opts.Pattern (1,1) string = "*_2D_*.dat"
        opts.Verbose (1,1) logical = true
    end

    files = dir(fullfile(folderPath, opts.Pattern));
    if isempty(files)
        error('No files match %s in "%s".', opts.Pattern, folderPath);
    end

    % Preallocate
    n = numel(files);
    Heatmap_data(n,1) = struct('name',"", 'H2O',NaN, 'H3O',NaN, ...
                               'temperature',NaN, 'sampling_freq',NaN, ...
                               'Position_1',"", 'Position_2',"", ...
                               'data',table());

    for i = 1:n
        filename = files(i).name;
        fp = fullfile(folderPath, filename);
        if opts.Verbose, fprintf("Loading: %s\n", filename); end

        % Load one file
        heatmap_tbl = io.import_heatmap_xyz(fp);

        % Parse name: <prefix>_<H2O>H2O_<H3O>H3O_<T>K_<fs>fs_2D_<AB>.dat
        % Example: AvOPs_130H2O_0H3O_300K_20fs_2D_XY.dat
        [~, raw_name] = fileparts(filename);
        tok = regexp(raw_name, '_(\d+)H2O_(\d+)H3O_(\d+)K_(\d+)fs_2D_([XYZ])([XYZ])$', 'tokens', 'once');
        
        if ~isempty(tok)
            H2O  = str2double(tok{1});
            H3O  = str2double(tok{2});
            T    = str2double(tok{3});
            fs   = str2double(tok{4});
            coord1 = string(tok{5});  % e.g., 'X'
            coord2 = string(tok{6});  % e.g., 'Y'
            
            struct_name = "Heatmap_" + tok{1} + "H2O_" + tok{2} + "H3O_" + ...
                         tok{3} + "K_" + tok{4} + "fs_2D_" + coord1 + coord2;
        else
            error("Could not parse filename format: %s\n" + ...
              	"Filename should be of format: \n" + ...
              	" - <prefix>_<H2O>H2O_<H3O>H3O_<T>K_<fs>fs_2D_<AB>.dat", filename);
        end

        % Fill struct
        Heatmap_data(i).name           = struct_name;
        Heatmap_data(i).H2O            = H2O;
        Heatmap_data(i).H3O            = H3O;
        Heatmap_data(i).temperature    = T;      % K
        Heatmap_data(i).sampling_freq  = fs;     % fs
        Heatmap_data(i).Position_1     = coord1; % 'X', 'Y', or 'Z'
        Heatmap_data(i).Position_2     = coord2; % 'X', 'Y', or 'Z'
        Heatmap_data(i).data           = heatmap_tbl;
    end
    
    T = struct2table(Heatmap_data);
    Mol_values = unique(T(:, {'H2O','H3O'}), 'rows');
end

