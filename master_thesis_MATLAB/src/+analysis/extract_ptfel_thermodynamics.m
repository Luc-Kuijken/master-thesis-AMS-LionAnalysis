function [PTFEL_params, PTFEL_data_processed] = extract_ptfel_thermodynamics(PTFEL_data, k_B, Av, opts)
%EXTRACT_PTFEL_THERMODYNAMICS Extract barriers and fit temperature dependence.  
%  [PTFEL_params, PTFEL_data_processed] = EXTRACT_PTFEL_THERMODYNAMICS(PTFEL_data, k_B, Av)
%
%  Extracts activation barriers from 1D PT-FEL curves and fits ΔH‡/ΔS‡.  
%  For water PT-FEL: single symmetric barrier
%  For MOF PT-FEL: two asymmetric barriers (left and right)
%
%  Inputs:
%    PTFEL_data : struct array from io.load_ptfel_struct (2fs only)
%    k_B        : Boltzmann constant [J/K]
%    Av         : Avogadro constant [1/mol]
%    opts.SavGolWindow : Savitzky-Golay window (default: 31)
%    opts.SavGolOrder  : Savitzky-Golay order (default: 3)
%
%  Outputs:
%    PTFEL_params : table with columns:
%                   loading, type, direction, dG_300K, dG_350K, dG_400K,
%                   dH, dH_kJmol, dS, dS_JmolK
%    PTFEL_data_processed : original struct with added barrier fields

    arguments
        PTFEL_data (1,:) struct
        k_B (1,1) double
        Av (1,1) double
        opts.SavGolWindow (1,1) double = 31
        opts.SavGolOrder (1,1) double = 3
    end
    
    % Filter for 1D data only
    mask_1d = ~[PTFEL_data.is_2d];
    PTFEL_1d = PTFEL_data(mask_1d);
    
    % Get unique loadings
    Table = struct2table(PTFEL_1d);
    unique_loadings = unique(Table.H2O + Table.H3O);
    
    % Preallocate results table
    n_loadings = numel(unique_loadings);
    result_rows = [];
    
    % Copy original data to processed version
    PTFEL_data_processed = PTFEL_data;
    
    % Process each loading
    for i = 1:n_loadings
        loading = unique_loadings(i);
        
        % Process water PT-FEL (symmetric, single barrier)
        water_params = process_water_ptfel(PTFEL_1d, loading, k_B, Av, opts);
        if ~isempty(water_params)
            result_rows = [result_rows; water_params];
        end
        
        % Process MOF PT-FEL (asymmetric, two barriers)
        mof_params = process_mof_ptfel(PTFEL_1d, loading, k_B, Av, opts);
        if ~isempty(mof_params)
            result_rows = [result_rows; mof_params];
        end
    end
    
    % Convert to table
    PTFEL_params = struct2table(result_rows);
    
    % Extract and store smoothed data for plotting
    for k = 1:numel(PTFEL_data_processed)
        % Only process 1D PT-FEL data
        if PTFEL_data_processed(k).is_2d
            continue
        end
        
        delta = PTFEL_data_processed(k).data.delta;
        fe = PTFEL_data_processed(k).data.free_energy;
        fe = fe - min(fe);  % Shift to zero minimum
        
        % Extract smoothed curve
        result = extract_ptfel_barrier_local(delta, fe, ...
                                            opts.SavGolWindow, ...
                                            opts.SavGolOrder);
        
        % Store smoothed data and spline function
        PTFEL_data_processed(k).fe_smooth = result.fe_smooth;
        PTFEL_data_processed(k).spline_func = result.spline_func;
        PTFEL_data_processed(k).barrier_kT = result.barrier_kT;
        PTFEL_data_processed(k).delta_saddle = result.delta_saddle;
        PTFEL_data_processed(k).fe_saddle = result.fe_saddle;
    end
end

%% Helper: Process water PTFEL (symmetric)
function params = process_water_ptfel(PTFEL_1d, loading, k_B, Av, opts)
    % Find all water PT-FEL for this loading
    mask = ([PTFEL_1d.type] == "water") & ...
           (([PTFEL_1d.H2O] + [PTFEL_1d.H3O]) == loading);
    
    if sum(mask) < 3
        params = [];
        return
    end
    
    water_data = PTFEL_1d(mask);
    temps = [water_data.temperature];
    
    % Extract barriers at each temperature
    dG_kT = zeros(1, numel(water_data));
    for k = 1:numel(water_data)
        delta = water_data(k).data.delta;
        fe = water_data(k).data.free_energy - min(water_data(k).data.free_energy);
        
        % Extract symmetric barrier
        result = extract_ptfel_barrier_local(delta, fe, ... 
                                            opts.SavGolWindow, ... 
                                            opts.SavGolOrder);
        dG_kT(k) = result.barrier_kT;
    end
    
    % Convert to Joules
    dG_J = dG_kT .* k_B .* temps;
    
    % Fit:    ΔG = ΔH - T·ΔS
    p = polyfit(temps, dG_J, 1);
    dS = -p(1);         % Slope = -ΔS [J/K]
    dH = p(2);          % Intercept = ΔH [J]
    dH_kJmol = dH * Av / 1000;
    dS_JmolK = dS * Av;
    
    % Store barriers at each temperature
    dG_300K = dH - 300 * dS;
    dG_350K = dH - 350 * dS;
    dG_400K = dH - 400 * dS;
    
    % Return struct
    params = struct('loading', loading, ... 
                    'type', "water", ... 
                    'direction', "symmetric", ...
                    'dG_300K_J', dG_300K, ...  
                    'dG_350K_J', dG_350K, ... 
                    'dG_400K_J', dG_400K, ...
                    'dG_300K_kJmol', dG_300K * Av / 1000, ... 
                    'dG_350K_kJmol', dG_350K * Av / 1000, ...
                    'dG_400K_kJmol', dG_400K * Av / 1000, ...
                    'dH_J', dH, ... 
                    'dH_kJmol', dH_kJmol, ...
                    'dS_JK', dS, ...
                    'dS_JmolK', dS_JmolK);
end

%% Helper:  Process MOF PT-FEL (asymmetric, two directions)
function params = process_mof_ptfel(PTFEL_1d, loading, k_B, Av, opts)
    % Find all MOF PT-FEL for this loading
    mask = ([PTFEL_1d.type] == "mof") & ...
           (([PTFEL_1d.H2O] + [PTFEL_1d.H3O]) == loading);
    
    if sum(mask) < 3
        params = [];
        return
    end
    
    mof_data = PTFEL_1d(mask);
    temps = [mof_data.temperature];
    
    % Extract left and right barriers at each temperature
    dG_left_kT = zeros(1, numel(mof_data));
    dG_right_kT = zeros(1, numel(mof_data));
    
    for k = 1:numel(mof_data)
        delta = mof_data(k).data.delta;
        fe = mof_data(k).data.free_energy - min(mof_data(k).data.free_energy);
        
        % Extract asymmetric barriers
        result = extract_ptfel_barrier_local(delta, fe, ...
                                            opts.SavGolWindow, ...
                                            opts.SavGolOrder);
        
        % Left barrier:    from left minimum to saddle
        dG_left_kT(k) = result.fe_saddle - result.fe_min_left;
        
        % Right barrier:  from right minimum to saddle
        dG_right_kT(k) = result.fe_saddle - result.fe_min_right;
    end
    
    % Fit left direction
    dG_left_J = dG_left_kT .* k_B .* temps;
    p_left = polyfit(temps, dG_left_J, 1);
    dS_left = -p_left(1);
    dH_left = p_left(2);
    
    % Fit right direction
    dG_right_J = dG_right_kT .* k_B .* temps;
    p_right = polyfit(temps, dG_right_J, 1);
    dS_right = -p_right(1);
    dH_right = p_right(2);
    
    % Return two structs (left and right)
    params_left = struct('loading', loading, ...
                         'type', "mof", ... 
                         'direction', "left", ...
                         'dG_300K_J', dH_left - 300 * dS_left, ...
                         'dG_350K_J', dH_left - 350 * dS_left, ...
                         'dG_400K_J', dH_left - 400 * dS_left, ...  
                         'dG_300K_kJmol', (dH_left - 300 * dS_left) * Av / 1000, ... 
                         'dG_350K_kJmol', (dH_left - 350 * dS_left) * Av / 1000, ... 
                         'dG_400K_kJmol', (dH_left - 400 * dS_left) * Av / 1000, ... 
                         'dH_J', dH_left, ...
                         'dH_kJmol', dH_left * Av / 1000, ... 
                         'dS_JK', dS_left, ... 
                         'dS_JmolK', dS_left * Av);
    
    params_right = struct('loading', loading, ...
                          'type', "mof", ...
                          'direction', "right", ...  
                          'dG_300K_J', dH_right - 300 * dS_right, ...
                          'dG_350K_J', dH_right - 350 * dS_right, ... 
                          'dG_400K_J', dH_right - 400 * dS_right, ... 
                          'dG_300K_kJmol', (dH_right - 300 * dS_right) * Av / 1000, ... 
                          'dG_350K_kJmol', (dH_right - 350 * dS_right) * Av / 1000, ...
                          'dG_400K_kJmol', (dH_right - 400 * dS_right) * Av / 1000, ...
                          'dH_J', dH_right, ... 
                          'dH_kJmol', dH_right * Av / 1000, ...
                          'dS_JK', dS_right, ...
                          'dS_JmolK', dS_right * Av);
    
    params = [params_left; params_right];
end

%% Helper: Extract PT-FEL barrier using Savitzky-Golay smoothing
function result = extract_ptfel_barrier_local(delta, free_energy, savgol_window, savgol_order)
%EXTRACT_PTFEL_BARRIER_LOCAL Extract activation barrier from 1D PT-FEL.

    % Convert to column vectors
    delta = delta(:);
    fe = free_energy(:);
    
    % Sort by delta (critical for Savitzky-Golay)
    [delta, idx] = sort(delta);
    fe = fe(idx);
    
    % Find indices where data is finite
    is_fin = isfinite(fe);
    fe_smooth = nan(size(fe)); % Preallocate with NaN
    
    if sum(is_fin) >= 5 % Minimum points required for a meaningful fit
        fe_finite = fe(is_fin);
        
        % Adjust window size to segment length
        window = min(savgol_window, floor(numel(fe_finite) / 2) * 2 - 1);
        if window < 5, window = 5; end
        if mod(window, 2) == 0, window = window - 1; end
        
        % Apply filter ONLY to finite segment
        smoothed_segment = sgolayfilt(fe_finite, savgol_order, window);
        
        % Map back to original indices to keep array lengths identical
        fe_smooth(is_fin) = smoothed_segment;
    else
        fe_smooth = fe; % Fallback if data is mostly Inf
    end
    
    % Create interpolation function (using only finite smoothed data)
    finite_mask = isfinite(fe_smooth);
    if sum(finite_mask) > 1
        spline_func = griddedInterpolant(delta(finite_mask), fe_smooth(finite_mask), 'spline', 'linear');
    else
        spline_func = [];
    end
    
    % Find saddle: maximum in central region (-0.5 to 0.5 Å)
    center_mask = (delta >= -0.5) & (delta <= 0.5) & is_fin;
    if sum(center_mask) < 2
        center_mask = is_fin; % Fallback to all finite points
    end
    
    [fe_saddle, saddle_idx_local] = max(fe_smooth(center_mask));
    delta_center = delta(center_mask);
    delta_saddle = delta_center(saddle_idx_local);
    
    % Find minima on left and right of saddle within finite range
    left_mask = (delta < delta_saddle) & is_fin;
    right_mask = (delta > delta_saddle) & is_fin;
    
    if sum(left_mask) < 1 || sum(right_mask) < 1
        barrier_kT = NaN;
        fe_min_left = NaN; fe_min_right = NaN;
        delta_min_left = NaN; delta_min_right = NaN;
    else
        [fe_min_left, left_idx_local] = min(fe_smooth(left_mask));
        delta_left = delta(left_mask);
        delta_min_left = delta_left(left_idx_local);
        
        [fe_min_right, right_idx_local] = min(fe_smooth(right_mask));
        delta_right = delta(right_mask);
        delta_min_right = delta_right(right_idx_local);
        
        % Barrier = saddle height - average of minima
        fe_min_avg = (fe_min_left + fe_min_right) / 2;
        barrier_kT = fe_saddle - fe_min_avg;
    end
    
    % Package results
    result = struct(...  
        'barrier_kT', barrier_kT, ... 
        'delta_saddle', delta_saddle, ...
        'fe_saddle', fe_saddle, ...
        'delta_min_left', delta_min_left, ...
        'fe_min_left', fe_min_left, ...
        'delta_min_right', delta_min_right, ...
        'fe_min_right', fe_min_right, ...
        'fe_smooth', fe_smooth, ... % Length matches delta
        'spline_func', spline_func);
end