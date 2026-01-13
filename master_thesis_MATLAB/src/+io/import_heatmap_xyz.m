function heatmap_data = import_heatmap_xyz(filename, dataLines)
%IMPORT_HEATMAP_XYZ Import 2D heatmap data from a text file
%  HEATMAP_DATA = IMPORT_HEATMAP_XYZ(FILENAME) reads data from text file
%  FILENAME for the default selection. Returns the data as a table.
%
%  HEATMAP_DATA = IMPORT_HEATMAP_XYZ(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  heatmap_data = import_heatmap_xyz("data/processed/heatmap_data/AvOPs_130H2O_0H3O_300K_20fs_2D_XY.dat", [3, Inf]);
%
%  See also READTABLE.

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [3, Inf];  % Skip first 2 header lines
end

%% Extract coordinate labels from filename
% Expected format: *_2D_AB.dat where A and B are coordinate labels (X, Y, Z)
[~, fname, ~] = fileparts(filename);
tok = regexp(fname, '_2D_([XYZ])([XYZ])$', 'tokens', 'once');

if ~isempty(tok)
    coord1 = tok{1};  % First coordinate (e.g., 'X')
    coord2 = tok{2};  % Second coordinate (e.g., 'Y')
    var1_name = "Position_" + lower(coord1);
    var2_name = "Position_" + lower(coord2);
else
    % Fallback if pattern not found
    warning('Could not parse coordinate labels from filename "%s". Using default names.', fname);
    var1_name = "Position_1";
    var2_name = "Position_2";
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = [var1_name, var2_name, "value"];
opts.VariableTypes = ["double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Import the data
heatmap_data = readtable(filename, opts);
% Remove empty rows
heatmap_data = heatmap_data(~all(ismissing(heatmap_data), 2), :);
end