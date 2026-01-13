function Diffusion_data = fit_arrhenius(Diffusion_data, Mol_values, k_B, Av)
%FIT_ARRHENIUS Fit Arrhenius params per H2O group; writes Ea [J], D0 [m^2/s].
%  Diffusion_data = FIT_ARRHENIUS(Diffusion_data, H2O_values, k_B)
%  Inputs: Diffusion_data(.H2O, .temperature[K], .D[m^2/s]); H2O_values; k_B[J/K]; Av[1/mol].
%  Output: Diffusion_data with .Ea [J], .Ea [kJ/mol] and .D0 [m^2/s] set for each group.
%  Requires ≥2 valid (T,D) pairs per group.

    for i = 1:numel(Mol_values.H2O)
        idx   = [Diffusion_data.H2O] == Mol_values.H2O(i);
        group = Diffusion_data(idx);

        T = [group.temperature];   % K
        D = [group.D];             % m^2/s

        % Need at least 2 points to fit
        if numel(T) < 2
            warning('Arrhenius fit skipped for %g H2O and %g H3O (need >=2 valid points).', Mol_values{i,:});
            [Diffusion_data(idx).Ea] = deal(NaN);
            [Diffusion_data(idx).D0] = deal(NaN);
            continue
        end

        % Sort by temperature
        [T, order] = sort(T(:));
        D = D(order);

        % Arrhenius fit: ln(D) = ln(D0) - (Ea/k_B)*(1/T)
        X = 1 ./ T;                 % 1/K
        Y = log(D);                 % ln(m^2/s)

        p  = polyfit(X, Y, 1);
        Ea = -p(1) * k_B;           % J
        D0 = exp(p(2));             % m^2/s

        % Store same group-level values in each member of this H2O group
        [Diffusion_data(idx).Ea] = deal(Ea);
        [Diffusion_data(idx).Ea_kJmol] = deal(Ea * Av / 1000);
        [Diffusion_data(idx).D0] = deal(D0);
    end
end

