function delta_l_loss = get_delta_l_loss(x)
    % Hutter: Equation (16)
    load('monopod_parameters', 'K_l')

    q = x(1:3);
    dq = x(4:6);
    dq_plus = get_impact(x);

    % Calculate the mass matrix
    M = get_M(q);

    % Calculate the system kinetic energy
    E_kin = 0.5 * dq' * M * dq;
    E_kin_plus = 0.5 * dq_plus' * M * dq_plus;
    
    % Calculate energy difference
    delta_E = E_kin_plus - E_kin;
    
    % Calculate delta_l_loss
    delta_l_loss = -sign(delta_E) * sqrt(2 / K_l * abs(delta_E));
end

