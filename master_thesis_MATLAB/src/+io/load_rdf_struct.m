function [RDF_data, Mol_values] = load_rdf_struct(folderPath, opts)
%LOAD_RDF_STRUCT Scan a folder of RDF files and build RDF_data struct. 
%  [RDF_data, Mol_values] = io.load_rdf_struct("data/rdf", Pattern="*.txt", Verbose=true)
%
%  Supports multiple RDF types: 
%    - rdf_o_h_*       : Water O-H (OH)
%    - rdf_o_o_*       : Water O-O (OO)
%    - rdf_linkero_h_* :  Linker O - Water H (LINKEROH)
%    - rdf_linkero_o_* :  Linker O - Water O (LINKEROO)
%    - rdf_mu3o_h_*    : Mu3-oxo O - Water H (MU3OH)
%    - rdf_mu3o_o_*    : Mu3-oxo O - Water O (MU3OO)
%
%  Expected filename format:
%    rdf_<type1>_<type2>_<H2O>H2O_<H3O>H3O_<T>K_<fs>fs.txt
%    Example: rdf_o_h_98H2O_0H3O_300K_100fs.txt
%
%  Outputs:
%    RDF_data   : struct array with fields:
%                   .name, .H2O, .H3O, .temperature, .sampling_freq, .rdf_type, .data
%    Mol_values :  table with unique H2O, H3O combinations

    arguments
        folderPath (1,1) string
        opts.Pattern (1,1) string = "rdf_*.txt"
        opts.Verbose (1,1) logical = true
    end

    files = dir(fullfile(folderPath, opts.Pattern));
    if isempty(files)
        error('No files match %s in "%s".', opts.Pattern, folderPath);
    end

    % Preallocate
    n = numel(files);
    RDF_data(n,1) = struct('name', "", 'H2O', NaN, 'H3O', NaN, ...
                           'temperature', NaN, 'sampling_freq', NaN, ...
                           'rdf_type', "", 'data', table());

    valid_count = 0;
    
    for i = 1:n
        filename = files(i).name;
        fp = fullfile(folderPath, filename);
        
        % Parse filename with flexible pattern for different RDF types
        [~, raw_name] = fileparts(filename);
        
        % Pattern:  rdf_<atom1>_<atom2>_<H2O>H2O_<H3O>H3O_<T>K_<fs>fs
        % Match various atom type combinations
        tok = regexp(raw_name, '^rdf_([a-z0-9]+)_([a-z])_(\d+)H2O_(\d+)H3O_(\d+)K_(\d+)fs$', 'tokens', 'once');
        
        if isempty(tok)
            if opts.Verbose
                warning("Could not parse filename: %s. Skipping.", filename);
            end
            continue
        end
        
        % Extract components
        atom1_str = tok{1};  % e.g., "o", "linkero", "mu3o"
        atom2_str = tok{2};  % e.g., "h", "o"
        H2O = str2double(tok{3});
        H3O = str2double(tok{4});
        T = str2double(tok{5});
        fs = str2double(tok{6});
        
        % Determine RDF type based on atom identifiers
        rdf_type = determine_rdf_type(atom1_str, atom2_str);
        
        if rdf_type == ""
            if opts.Verbose
                warning("Unknown RDF type for: %s. Skipping.", filename);
            end
            continue
        end
        
        if opts.Verbose
            fprintf("Loading:  %s (%s, %dH2O, %dH3O, %dK)\n", filename, rdf_type, H2O, H3O, T);
        end
        
        % Load RDF data
        rdf_tbl = io.import_rdf(fp);
        
        % Fill struct
        valid_count = valid_count + 1;
        struct_name = sprintf("rdf_%s_%dH2O_%dH3O_%dK_%dfs", rdf_type, H2O, H3O, T, fs);
        
        RDF_data(valid_count).name = struct_name;
        RDF_data(valid_count).H2O = H2O;
        RDF_data(valid_count).H3O = H3O;
        RDF_data(valid_count).temperature = T;
        RDF_data(valid_count).sampling_freq = fs;
        RDF_data(valid_count).rdf_type = rdf_type;
        RDF_data(valid_count).data = rdf_tbl;
    end
    
    % Trim to valid entries
    RDF_data = RDF_data(1:valid_count);
    
    % Get unique molecule values
    T = struct2table(RDF_data);
    Mol_values = unique(T(: , {'H2O', 'H3O'}), 'rows');
    Mol_values.Total = Mol_values.H2O + Mol_values.H3O;
end

%% Helper function to determine RDF type
function rdf_type = determine_rdf_type(atom1_str, atom2_str)
    % Normalize to uppercase for comparison
    atom1 = upper(atom1_str);
    atom2 = upper(atom2_str);
    
    % Map atom combinations to RDF types
    if strcmp(atom1, "O") && strcmp(atom2, "H")
        rdf_type = "O-H";
    elseif strcmp(atom1, "O") && strcmp(atom2, "O")
        rdf_type = "O-O";
    elseif strcmp(atom1, "LINKERO") && strcmp(atom2, "H")
        rdf_type = "LinkerO-H";
    elseif strcmp(atom1, "LINKERO") && strcmp(atom2, "O")
        rdf_type = "LinkerO-O";
    elseif strcmp(atom1, "MU3O") && strcmp(atom2, "H")
        rdf_type = "Mu3O-H";
    elseif strcmp(atom1, "MU3O") && strcmp(atom2, "O")
        rdf_type = "Mu3O-O";
    else
        rdf_type = "";  % Unknown type
    end
end