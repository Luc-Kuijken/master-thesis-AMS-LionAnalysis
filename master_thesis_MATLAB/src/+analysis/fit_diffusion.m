function Diffusion_data = fit_diffusion(Diffusion_data, start_time, end_time)
%FIT_DIFFUSION Estimate D from MSD=6*D*t+c over [start_time, end_time].
%  Diffusion_data = FIT_DIFFUSION(Diffusion_data, start_time, end_time)
%  Inputs: 
%    Diffusion_data(.data.time[ps], .data.msd_a[Angstrom^2])
%    start_time: starting time in ps
%    end_time: ending time in ps
%  Output: adds .D [m^2/s], .fit_coeffs [slope,intercept], .fit_window [i1 i2].

    for k = 1:numel(Diffusion_data)
        time_ps = Diffusion_data(k).data.time;     % [ps]
        msd_a2  = Diffusion_data(k).data.msd_a;    % [Angstrom^2]

        N = numel(time_ps);
        
        % Find closest indices to start_time and end_time
        [~, i1] = min(abs(time_ps - start_time));
        [~, i2] = min(abs(time_ps - end_time));
        
        % Ensure valid range
        i1 = max(1, min(N-1, i1));
        i2 = max(i1+1, min(N, i2));

        tSlice = time_ps(i1:i2);
        ySlice = msd_a2(i1:i2);

        % Fit the linear diffusion line
        p = polyfit(tSlice, ySlice, 1);            % slope [Angstrom^2/ps]
        D = (p(1) / 6) * (1e-20 / 1e-12);          % -> [m^2/s]

        Diffusion_data(k).D           = D;
        Diffusion_data(k).fit_coeffs  = p;
        Diffusion_data(k).fit_window  = [i1 i2];
    end
end

