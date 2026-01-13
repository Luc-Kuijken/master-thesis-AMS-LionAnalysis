function [PTFEL_data, Mol_values] = load_ptfel_struct(folderPath, opts)
%LOAD_PTFEL_STRUCT Scan folder for PT-FEL files and build PTFEL_data struct.
%  [PTFEL_data, Mol_values] = io.load_ptfel_struct("data/ptfel", Pattern="*.txt", Verbose=true)
%
%  Expected filename format:
%    ptfel_water_H2O_H3O_K_fs.txt      (1D, 3 columns)
%    ptfel_water_2d_H2O_H3O_K_fs.txt   (2D, 4 columns)
%    ptfel_mof_H2O_H3O_K_fs.txt        (1D, 3 columns)
%    ptfel_mof_2d_H2O_H3O_K_fs.txt     (2D, 4 columns)
%
%  Outputs:
%    PTFEL_data   : struct array with fields:
%                   .name, .H2O, .H3O, .temperature, .sampling_freq, 
%                   .type ("water" or "mof"), .is_2d (logical), .data
%    Mol_values   : table with unique H2O, H3O combinations

    arguments
        folderPath (1,1) string
        opts.Pattern (1,1) string = "ptfel_*.txt"
        opts.SamplingFreq (1,: ) double = []  % Filter by sampling frequency (e.g., [2])
        opts.Verbose (1,1) logical = true
    end

    files = dir(fullfile(folderPath, opts.Pattern));
    if isempty(files)
        error('No files match %s in "%s".', opts.Pattern, folderPath);
    end

    % Preallocate
    n = numel(files);
    PTFEL_data(n,1) = struct('name', "", 'H2O', NaN, 'H3O', NaN, ... 
                            'temperature', NaN, 'sampling_freq', NaN, ...
                            'type', "", 'is_2d', false, 'data', table());

    valid_count = 0;
    
    for i = 1:n
        filename = files(i).name;
        fp = fullfile(folderPath, filename);
        
        % Parse filename:   ptfel_TYPE_H2O_H3O_K_fs.txt
        % or:              ptfel_TYPE_2d_H2O_H3O_K_fs.txt
        [~, raw_name] = fileparts(filename);
        
        % Try 2D pattern first
        tok = regexp(raw_name, '^ptfel_([a-z]+)_2d_(\d+)H2O_(\d+)H3O_(\d+)K_(\d+)fs$', 'tokens', 'once');
        
        if ~isempty(tok)
            % 2D file
            type_str = string(tok{1});  % "water" or "mof"
            H2O = str2double(tok{2});
            H3O = str2double(tok{3});
            T = str2double(tok{4});
            fs = str2double(tok{5});
            is_2d = true;
        else
            % Try 1D pattern
            tok = regexp(raw_name, '^ptfel_([a-z]+)_(\d+)H2O_(\d+)H3O_(\d+)K_(\d+)fs$', 'tokens', 'once');
            
            if ~isempty(tok)
                % 1D file
                type_str = string(tok{1});
                H2O = str2double(tok{2});
                H3O = str2double(tok{3});
                T = str2double(tok{4});
                fs = str2double(tok{5});
                is_2d = false;
            else
                if opts.Verbose
                    warning("Could not parse filename: %s. Skipping.", filename);
                end
                continue
            end
        end
        
        % Apply sampling frequency filter
        if ~isempty(opts.SamplingFreq) && ~ismember(fs, opts.SamplingFreq)
            if opts.Verbose
                fprintf("Skipping:  %s (sampling freq %dfs not in filter)\n", filename, fs);
            end
            continue
        end
        
        if opts.Verbose
            if is_2d
                fprintf("Loading:   %s (%s, 2D, %dH2O, %dH3O, %dK, %dfs)\n", ...
                       filename, type_str, H2O, H3O, T, fs);
            else
                fprintf("Loading:  %s (%s, 1D, %dH2O, %dH3O, %dK, %dfs)\n", ...
                       filename, type_str, H2O, H3O, T, fs);
            end
        end
        
        % Load PTFEL data (automatically handles 1D or 2D format)
        ptfel_tbl = io.import_ptfel(fp);
        
        % Fill struct
        valid_count = valid_count + 1;
        
        if is_2d
            struct_name = sprintf("ptfel_%s_2d_%dH2O_%dH3O_%dK_%dfs", ...
                                 type_str, H2O, H3O, T, fs);
        else
            struct_name = sprintf("ptfel_%s_%dH2O_%dH3O_%dK_%dfs", ...
                                 type_str, H2O, H3O, T, fs);
        end
        
        PTFEL_data(valid_count).name = struct_name;
        PTFEL_data(valid_count).H2O = H2O;
        PTFEL_data(valid_count).H3O = H3O;
        PTFEL_data(valid_count).temperature = T;
        PTFEL_data(valid_count).sampling_freq = fs;
        PTFEL_data(valid_count).type = type_str;
        PTFEL_data(valid_count).is_2d = is_2d;
        PTFEL_data(valid_count).data = ptfel_tbl;
    end
    
    % Trim to valid entries
    PTFEL_data = PTFEL_data(1:valid_count);
    
    % Get unique molecule values
    T = struct2table(PTFEL_data);
    Mol_values = unique(T(: , {'H2O', 'H3O'}), 'rows');
    Mol_values.Total = Mol_values.H2O + Mol_values.H3O;
end