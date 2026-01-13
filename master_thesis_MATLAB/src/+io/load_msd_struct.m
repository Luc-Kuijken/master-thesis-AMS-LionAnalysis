function [Diffusion_data, Mol_values] = load_msd_struct(folderPath, opts)
%LOAD_MSD_STRUCT Scan a folder of MSD files and build Diffusion_data.
%  [D, H2O_vals] = io.load_msd_folder("data/raw", Pattern="*.dat", Smooth=true)
%  Also supports *.txt MSD files with Pattern="*.txt"

    arguments
        folderPath (1,1) string
        opts.Pattern (1,1) string = "*.txt"
        opts.Smooth  (1,1) logical = true
        opts.SmoothWindow (1,1) double = NaN
        opts.Verbose (1,1) logical = true
    end

    files = dir(fullfile(folderPath, opts.Pattern));
    if isempty(files)
        error('No files match %s in "%s".', opts.Pattern, folderPath);
    end

    % Preallocate
    n = numel(files);
    Diffusion_data(n,1) = struct('name',"", 'H2O',NaN, 'H3O',NaN, ...
                                 'temperature',NaN, 'sampling',NaN, ...
                                 'data',table());

    for i = 1:n
        filename = files(i).name;
        fp = fullfile(folderPath, filename);
        if opts.Verbose, fprintf("Loading: %s\n", filename); end

        % Load one file based on file extension
        [~, ~, ext] = fileparts(filename);
        if strcmpi(ext, '.txt')
            msd = io.import_msd_txt(fp);
        elseif strcmpi(ext, '.dat')
            msd = io.import_msd(fp);
        else
            warning("Skipping %s... Unknown file extension " + ...
              	"(Only uses .txt and .dat files).", filename)
            continue
        end

        % Optional smoothing for derived signals
        if opts.Smooth
            if isnan(opts.SmoothWindow)
                win = min(1001, floor(numel(msd.msd_a)/100)+1);
            else
                win = max(1, round(opts.SmoothWindow));
                if mod(win,2)==0, win = win+1; end
            end
            avg_msd = movmean(msd.msd_a, win);
        else
            avg_msd = msd.msd_a;
            win = 11;
        end

        % Derived msd properties
        % dmsd/dt in m^2/s
        msd.dmsd_dt = (gradient(avg_msd) ./ gradient(msd.time)) * (1e-20/1e-12); % m^2/s
        
        % Beta: d(ln MSD)/d(ln t)
        log_msd = log(avg_msd);
        log_time = log(msd.time);
        % Compute gradient
        beta = gradient(log_msd) ./ gradient(log_time);
        msd.beta = movmean(beta, 51);

        % Parse name: <prefix>_<H2O>H2O_<H3O>H3O_<T>K.dat or msd_<H2O>H2O_<H3O>H3O_<T>K_*.txt
        [~, raw_name] = fileparts(filename);
        tok = regexp(raw_name, '_(\d+)H2O_(\d+)H3O_(\d+)K', 'tokens', 'once');
        tok2 = regexp(raw_name, '_(\d+)fs', 'tokens', 'once');
        
        if ~isempty(tok)
            H2O = str2double(tok{1}); 
            H3O = str2double(tok{2}); 
            T = str2double(tok{3});
            
            struct_name = "msd_" + tok{1} + "H2O_" + tok{2} + "H3O_" + tok{3} + "K";
            if ~isempty(tok2)
                fs = str2double(tok2{1}); 
            
                struct_name = struct_name + "_" + tok2{1} + "fs";
            else
                fs = NaN;
            end
        else
            error("Could not parse filename format: %s\n" + ...
              	"Filename should be one of either format: \n" + ...
              	" - <prefix>_<H2O>H2O_<H3O>H3O_<T>K.dat \n" + ...
              	" - <prefix>_<H2O>H2O_<H3O>H3O_<T>K_*.txt", filename);
        end

        % Fill struct
        Diffusion_data(i).name        = struct_name;
        Diffusion_data(i).H2O         = H2O;                                % Number of water molecules
        Diffusion_data(i).H3O         = H3O;                                % Number of hydronium molecules
        Diffusion_data(i).temperature = T;                                  % Temperature in K
        Diffusion_data(i).sampling    = fs;                                	% Sampling frequency in fs 
        Diffusion_data(i).data        = msd;
    end
    
    T = struct2table(Diffusion_data);
    Mol_values = unique(T(:, {'H2O','H3O'}), 'rows');
    Mol_values.Total = Mol_values.H2O + Mol_values.H3O;
end

