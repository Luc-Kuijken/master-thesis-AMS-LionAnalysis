function ptfel = import_ptfel(filename, dataLines)
%IMPORT_PTFEL Import PT-FEL data from a text file (1D or 2D format).
%  ptfel = IMPORT_PTFEL(FILENAME) reads PT-FEL data from text file.
%  
%  Automatically detects file format: 
%    - 1D format: 3 columns (delta, free_energy, count)
%    - 2D format: 4 columns (delta, distance, free_energy, count)
%
%  Returns table with appropriate column names.
%
%  Example:
%  ptfel = io.import_ptfel("data/ptfel/ptfel_water_126H2O_4H3O_400K_2fs.txt");

    arguments
        filename (1,1) string
        dataLines (1,2) double = [1, Inf]
    end
    
    % Read first valid data line to determine format
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    
    % Skip lines until we find first data line
    is_2d = false;
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && ~isempty(line) && line(1) ~= '#'
            % Split line and count columns
            parts = strsplit(strtrim(line));
            n_cols = numel(parts);
            
            if n_cols == 4
                is_2d = true;
            elseif n_cols == 3
                is_2d = false;
            else
                error('Unexpected number of columns (%d) in file: %s', n_cols, filename);
            end
            break
        end
    end
    fclose(fid);
    
    % Set up import options based on detected format
    if is_2d
        % 2D format: delta, distance, free_energy, count
        opts = delimitedTextImportOptions("NumVariables", 4);
        opts.VariableNames = ["delta", "distance", "free_energy", "count"];
        opts.VariableTypes = ["double", "double", "double", "double"];
    else
        % 1D format: delta, free_energy, count
        opts = delimitedTextImportOptions("NumVariables", 3);
        opts.VariableNames = ["delta", "free_energy", "count"];
        opts.VariableTypes = ["double", "double", "double"];
    end
    
    % Common import settings
    opts.DataLines = dataLines;
    opts.Delimiter = " ";
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts.ConsecutiveDelimitersRule = "join";
    opts.LeadingDelimitersRule = "ignore";
    opts.TrailingDelimitersRule = "ignore";
    opts.CommentStyle = "#";
    
    % Import the data
    ptfel = readtable(filename, opts);
    
    % Remove empty rows
    ptfel = ptfel(~all(ismissing(ptfel), 2), :);
end