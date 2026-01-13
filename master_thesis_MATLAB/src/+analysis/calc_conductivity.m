function Conductivity_data = calc_conductivity(Arrhenius_params, PTFEL_params, opts)
%CALC_CONDUCTIVITY Calculate vehicular and structural conductivity.   
%  Conductivity_data = CALC_CONDUCTIVITY(Arrhenius_params, PTFEL_params, opts)
%
%  Uses:    
%    sigma_veh(T) = (n q^2 D_0)/(k_B T) exp(-E_a / k_B T)
%    sigma_struct(T) = κ (n q^2 a^2)/(2dh) exp(-DeltaH/k_B T + DeltaS/k_B)
%
%  Inputs:    
%    Arrhenius_params : table with columns:  Total, Ea, Ea_kJmol, D0, T_avg, D_avg, RMSE
%    PTFEL_params     : table with dH_J, dS_JK, loading, type
%
%  Name-Value options:
%    V_cell        : Unit cell volume [Angstrom^3] (default: 9.0749e+03)
%    a_hop         : Hopping distance [Angstrom] (default: 2.45)
%    kappa         : Proton availability factor (default: 1.0)
%    d             : Dimensionality (default: 3)
%    temperatures  : Temperature range [K] (default: 300:10:400)
%    H3O_per_cell  : Number of H3O per unit cell (default: derived from loading)

    arguments
        Arrhenius_params (:,: ) table
        PTFEL_params (:,13) table
        opts.V_cell (1,1) double = 9.0749e+03
        opts.a_hop (1,1) double = 2.45
        opts.kappa (1,1) double = 1.0
        opts.d (1,1) double = 3
        opts.k_B (1,1) double = 1.380649e-23
        opts.h (1,1) double = 6.62607015e-34
        opts.e (1,1) double = 1.602176634e-19
        opts.temperatures (1,: ) double = 300:10:400
        opts.H3O_per_cell (1,: ) double = []  % Optional: specify H3O count per loading
    end
    
    % Convert units
    V_cell_m3 = opts.V_cell * 1e-30;
    a_hop_m = opts.a_hop * 1e-10;
    
    % Extract loadings from Arrhenius_params (uses 'Total' column)
    if ~ismember('Total', Arrhenius_params.Properties.VariableNames)
        error('Arrhenius_params must contain a "Total" column (H2O + H3O)');
    end
    
    loadings = Arrhenius_params.Total;
    n_loadings = numel(loadings);
    n_temps = numel(opts.temperatures);
    n_rows = n_loadings * n_temps;
    
    % Assume H3O count = 4
    if isempty(opts.H3O_per_cell)
        H3O_counts = 4;
        warning('H3O_per_cell not specified.Using default: 4 H2O');
    else
        if numel(opts.H3O_per_cell) == 1
            H3O_counts = repmat(opts.H3O_per_cell, n_loadings, 1);
        elseif numel(opts.H3O_per_cell) == n_loadings
            H3O_counts = opts.H3O_per_cell(: );
        else
            error('H3O_per_cell must be scalar or have same length as number of loadings');
        end
    end
    
    % Preallocate output table
    Loading = zeros(n_rows, 1);
    Temperature = zeros(n_rows, 1);
    sigma_veh = zeros(n_rows, 1);
    sigma_struct = zeros(n_rows, 1);
    sigma_total = zeros(n_rows, 1);
    
    sigma0_veh = zeros(n_rows, 1);
    sigma0_struct = zeros(n_rows, 1);
    dH_w = zeros(n_rows, 1);
    dS_w = zeros(n_rows, 1);
    
    row_idx = 1;
    
    for i = 1:n_loadings
        loading = loadings(i);
        H3O = H3O_counts(i);
        
        % Get Arrhenius parameters for this loading
        D0 = Arrhenius_params.D0(i);
        Ea = Arrhenius_params.Ea(i);
        
        % Number density of protons
        n_protons = H3O / V_cell_m3;
        
        % Get PTFEL parameters for water-water hopping
        water_mask = (PTFEL_params.loading == loading) & ...
                    (PTFEL_params.type == "water") & ...
                    (PTFEL_params.direction == "symmetric");
        
        if any(water_mask)
            dH_water = PTFEL_params.dH_J(water_mask);
            dS_water = PTFEL_params.dS_JK(water_mask);
        else
            dH_water = NaN;
            dS_water = NaN;
        end
        
        % Calculate conductivities at each temperature
        for j = 1:n_temps
            T = opts.temperatures(j);
            
            % Vehicular conductivity:
            sigma0_v = (n_protons * opts.e^2 * D0) / (opts.k_B);
            sigma_v = sigma0_v / T * exp(-Ea / (opts.k_B * T));
            
            % Structural conductivity (water-water):
            sigma0_s = opts.kappa * (n_protons * opts.e^2 * a_hop_m^2) / ...
                (2 * opts.d * opts.h);
            
            if ~isnan(dH_water)
                sigma_sw = sigma0_s * exp(-dH_water / (opts.k_B * T) + dS_water / opts.k_B);
            else
                sigma_sw = NaN;
            end
            
            % Total conductivity:
            sigma_t = sigma_v + (isfinite(sigma_sw) * sigma_sw);
            
            % Store in arrays
            Loading(row_idx) = loading;
            Temperature(row_idx) = T;
            sigma_veh(row_idx) = sigma_v;
            sigma_struct(row_idx) = sigma_sw;
            sigma_total(row_idx) = sigma_t;
            
            sigma0_veh(row_idx) = sigma0_v;
            sigma0_struct(row_idx) = sigma0_s;
            dH_w(row_idx) = dH_water;
            dS_w(row_idx) = dS_water;
            
            row_idx = row_idx + 1;
        end
    end
    
    % Convert to table
    Conductivity_data = table(Loading, Temperature, sigma_veh, ...
        sigma_struct, sigma_total, ...
        sigma0_veh, sigma0_struct, dH_w, dS_w);
end