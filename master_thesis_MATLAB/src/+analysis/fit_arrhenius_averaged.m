function Arrhenius_params = fit_arrhenius_averaged(Diffusion_data, Mol_values, k_B, Av)
%FIT_ARRHENIUS_AVERAGED Fit Arrhenius parameters per loading by averaging D at each T.
%  Arrhenius_avg = FIT_ARRHENIUS_AVERAGED(Diffusion_data, Mol_values, k_B, Av)
%  Inputs: 
%    Diffusion_data(.H2O,.H3O,.temperature[K],.D[m^2/s])
%    Mol_values table with columns H2O, H3O, Total
%    k_B [J/K]; Av [1/mol]
%  Output: 
%    Arrhenius_params table with columns:
%      Total, Ea [J], Ea_kJmol [kJ/mol], D0 [m^2/s], T_avg [K], D_avg [m^2/s],
%      RMSE [dimensionless in ln(D)]

    % Get unique loadings
    unique_loadings = unique(Mol_values.Total);
    n_loadings = numel(unique_loadings);
    
    % Preallocate output table
    Arrhenius_params = table('Size', [n_loadings, 7], ...
    	'VariableTypes', ...
        {'double', 'double', 'double',   'double', 'cell',  'cell',  'double'}, ...
    	'VariableNames', ...
        {'Total',  'Ea',     'Ea_kJmol', 'D0',     'T_avg', 'D_avg', 'RMSE'});
    
    for i = 1:n_loadings
        loading_val = unique_loadings(i);
        
        % Find all data with this loading
        idx = arrayfun(@(x) x.H2O + x.H3O, Diffusion_data) == loading_val;
        group = Diffusion_data(idx);
        
        % Get unique temperatures for this loading
        temps = unique([group.temperature]);
        n_temps = numel(temps);
        
        % Preallocate
        T = zeros(n_temps, 1);
        D_avg = zeros(n_temps, 1);
        
        % Average D at each temperature
        for ti = 1:n_temps
            temp_val = temps(ti);
            temp_mask = [group.temperature] == temp_val;
            D_values = [group(temp_mask).D];
            
            T(ti) = temp_val;
            D_avg(ti) = mean(D_values);
        end
        
        % Need at least 2 points to fit
        if numel(T) < 2
            warning('Not enough temperatures for loading %g (need >=2). Skipping.', loading_val);
            Arrhenius_params.Total(i) = loading_val;
            Arrhenius_params.Ea(i) = NaN;
            Arrhenius_params.Ea_kJmol(i) = NaN;
            Arrhenius_params.D0(i) = NaN;
            Arrhenius_params.T_avg{i} = [];
            Arrhenius_params.D_avg{i} = [];
            Arrhenius_params.RMSE(i) = NaN;
            continue
        end
        
        % Sort by temperature
        [T, order] = sort(T);
        D_avg = D_avg(order);
        
        % Arrhenius fit: ln(D) = ln(D0) - (Ea/k_B)*(1/T)
        X = 1 ./ T(:);          % 1/K
        Y = log(D_avg(:));      % ln(m^2/s)
        
        [p, S] = polyfit(X, Y, 1);
        Ea = -p(1) * k_B;           % J
        D0 = exp(p(2));             % m^2/s
        
        % Compute RMSE of fit (standard deviation of residuals)
        dof = S.df;
        RMSE = S.normr / sqrt(dof);
        
        % Store in table (including averaged points)
        Arrhenius_params.Total(i) = loading_val;
        Arrhenius_params.Ea(i) = Ea;
        Arrhenius_params.Ea_kJmol(i) = Ea * Av / 1000;
        Arrhenius_params.D0(i) = D0;
        Arrhenius_params.T_avg{i} = T;
        Arrhenius_params.D_avg{i} = D_avg;
        Arrhenius_params.RMSE(i) = RMSE;
    end
end