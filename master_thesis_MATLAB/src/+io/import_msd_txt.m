function msd = import_msd_txt(filename, dataLines)
%IMPORT_MSD_TXT Import data from a text file with MSD data in TXT format
%  MSD = IMPORT_MSD_TXT(FILENAME) reads data from text file
%  FILENAME for the default selection. Returns the data as a table.
%
%  MSD = IMPORT_MSD_TXT(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  msd = import_msd_txt("data/raw/msd_130H2O_0H3O_400K_100fs_water.txt");
%
%  See also READTABLE.

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    % Skip header lines (lines starting with '#')
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    
    lineCount = 0;
    headerLines = 0;
    while ~feof(fid)
        line = fgetl(fid);
        lineCount = lineCount + 1;
        if ~isempty(line) && (line(1) == '#')
            headerLines = lineCount;
        else
            break;
        end
    end
    fclose(fid);
    
    dataLines = [headerLines + 1, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["time", "msd_a", "sum", "count"];
opts.SelectedVariableNames = ["time", "msd_a"];
opts.VariableTypes = ["double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";
opts.TrailingDelimitersRule = "ignore";
opts.CommentStyle = "#";

% Import the data
msd = readtable(filename, opts);

end