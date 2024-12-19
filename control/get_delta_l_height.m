function delta_l_height = get_delta_l_height(x, H_d)
    % Derived from Raibert
    load('monopod_parameters', 'K_l')
    
    % Total system energy at low point
    E_low = get_E_sys(x);

    % Total system energy at high point
    q = x(1:3);
    E_high = get_E_sys_d(q, H_d);

    % Calculate delta_l_high
    delta_E = E_high - E_low;
    delta_l_height = sign(delta_E) * sqrt(2 / K_l * abs(delta_E));
end

