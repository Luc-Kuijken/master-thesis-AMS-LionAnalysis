function rdf = import_rdf(filename, dataLines)
%IMPORT_RDF Import RDF data from a text file.
%  rdf = IMPORT_RDF(FILENAME) reads RDF data from FILENAME.
%  Returns table with columns: r, rdf, irdf, pmf, rawcount, rawcount_norm
%
%  Expected format (space-delimited, # comments):
%    #r RDF iRDF PMF(kTunits) rawcount rawcountpertimestepperatom

    arguments
        filename (1,1) string
        dataLines (1,2) double = [2, Inf]  % Skip header line
    end

    % Set up import options
    opts = delimitedTextImportOptions("NumVariables", 6);
    
    opts.DataLines = dataLines;
    opts.Delimiter = " ";
    
    opts.VariableNames = ["r", "rdf", "irdf", "pmf", "rawcount", "rawcount_norm"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];
    
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts.ConsecutiveDelimitersRule = "join";
    opts.LeadingDelimitersRule = "ignore";
    opts.CommentStyle = "#";
    
    % Import data
    rdf = readtable(filename, opts);
    
    % Remove empty/NaN rows
    rdf = rdf(~all(ismissing(rdf), 2), :);
end