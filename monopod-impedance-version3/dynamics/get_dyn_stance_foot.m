function dx = get_dyn_stance_foot(x, l0_s, phi0_s, m, I1, I2, g, K_l, D_l, K_phi, D_phi)
    import casadi.*;
    dx = MX.zeros(6,1);

    % Compute the current leg length and hip angle
    r_F = get_foot_pos(x, 0.5, 0.5);
    l = norm(r_F);
    phi = atan2(r_F(2), r_F(1));

    % Compute impedance forces (spring-damper system)
    F_leg = K_l * (l0_s - l) - D_l * x(4);
    tau_hip = K_phi * (phi0_s - phi) - D_phi * x(5);

    % Compute joint torques using the Jacobian
    J_pq = get_J_pq(x);
    tau_imp = J_pq' * [F_leg; tau_hip];

    % Compute stance force (Ground Reaction Force)
    k_leg = 1000; % stiffness
    d_leg = 10;   % Added damping
    F_stance = -k_leg * (x(1) - l0_s) - d_leg * x(4);

    % Compute system dynamics
    dx(1) = x(4); % Vertical velocity
    dx(2) = x(5); % Joint 1 velocity
    dx(3) = x(6); % Joint 2 velocity
    dx(4) = (F_stance / m) - g; % Vertical acceleration (including ground reaction)
    dx(5) = tau_imp(1) / I1; % Joint 1 angular acceleration
    dx(6) = tau_imp(2) / I2; % Joint 2 angular acceleration
end
