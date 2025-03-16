function dx = get_dyn_stance_imp(x, l0, phi0, K_l, D_l, K_phi, D_phi, m, I1, I2, g, L1, L2)
    import casadi.*

    dx = MX.zeros(6,1);

    dx(1) = x(4);
    dx(2) = x(5);
    dx(3) = x(6);

    % Possibly a "stance spring" for body y. 
    % That’s optional: you can keep or remove it.
    k_leg = 2000;  % example
    F_leg_body = -k_leg*( x(1) - 0 ); % old style "foot on ground" 
    dd_y_body = F_leg_body/m - g;

    % Now also do polar leg impedance for the joints:
    [l, phi, dl, dphi] = compute_l_phi_and_derivs(x, L1, L2);

    F_l    = K_l   * (l0   - l)   - D_l   * dl;
    F_phi  = K_phi * (phi0 - phi) - D_phi * dphi;

    J_p      = get_polar_jacobian(x, L1, L2);
    tau_imp  = J_p.' * [F_l; F_phi];

    tau_hip  = tau_imp(1);
    tau_knee = tau_imp(2);

    dd_q1 = tau_hip  / I1;
    dd_q2 = tau_knee / I2;

    dx(4) = dd_y_body;
    dx(5) = dd_q1;
    dx(6) = dd_q2;
end
