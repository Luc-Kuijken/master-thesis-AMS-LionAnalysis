function msd = import_msd(filename, dataLines)
%IMPORT_MSD Import data from a text file
%  MSD = IMPORT_MSD(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  MSD = IMPORT_MSD(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  msd = import_msd("C:\Users\20192602\MATLAB Drive\Year 6\Master Thesis\msd_data\msd_130H2O_0H3O_400K.dat", [2, Inf]);
%
%  See also READTABLE.

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["time", "msd_a", "Var3"];
opts.SelectedVariableNames = ["time", "msd_a"];
opts.VariableTypes = ["double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, "Var3", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var3", "EmptyFieldRule", "auto");

% Import the data
msd = readtable(filename, opts);

end