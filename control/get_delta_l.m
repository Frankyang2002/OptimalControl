function delta_l_height = get_delta_l(x_prev, H_d)
    load('monopod_parameters', 'K_l', 'K_p')
    
    % Previous gravitational potential energy at high point
    E_prev = get_E_sys(x_prev);

    % Desired total system energy at high point
    q_prev = x_prev(1:3);
    E_d = get_E_sys_d(q_prev, H_d);

    % Calculate energy difference
    delta_E =  K_p * (E_d - E_prev);

    % Calculate delta_l_high    
    delta_l_height = sign(delta_E) * sqrt(2 / K_l * abs(delta_E));
end

