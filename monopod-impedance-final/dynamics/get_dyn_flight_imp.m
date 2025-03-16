function dx = get_dyn_flight_imp(x, l0, phi0, K_l, D_l, K_phi, D_phi, m, I1, I2, g, L1, L2)
    % x = [y; q1; q2; dy; dq1; dq2]
    % l0, phi0 are the "rest" length & angle that we are optimizing
    % K_l, D_l, etc. are the impedance gains you choose

    import casadi.*

    dx = MX.zeros(6,1);

    % 1) Positions
    dx(1) = x(4);   % y-dot
    dx(2) = x(5);   % q1-dot
    dx(3) = x(6);   % q2-dot

    % 2) Evaluate the polar coordinates (l, phi) and their time derivatives (dl, dphi)
    [l, phi, dl, dphi] = compute_l_phi_and_derivs(x, L1, L2);

    % 3) Impedance Force in polar coords
    F_l    = K_l   * (l0   - l)   - D_l   * dl;   % radial
    F_phi  = K_phi * (phi0 - phi) - D_phi * dphi; % angular

    % 4) Map to joint torques
    J_p        = get_polar_jacobian(x, L1, L2);  % 2x2, partial(r,phi)/partial(q1,q2)
    tau_imp    = J_p.' * [F_l; F_phi];

    tau_hip    = tau_imp(1); % joint1 torque
    tau_knee   = tau_imp(2); % joint2 torque

    % 5) Acceleration
    dd_y   = -g;   % no ground contact => only gravity
    dd_q1  = tau_hip  / I1;
    dd_q2  = tau_knee / I2;

    % 6) Return dx
    dx(4) = dd_y;
    dx(5) = dd_q1;
    dx(6) = dd_q2;
end
