function u = impedance_controller(x, x_set, control_mode)
    load('monopod_parameters.mat', 'K_l', 'K_phi', 'D_l', 'D_phi', 'l_UL', 'l_LL', 'enable_trq_limit', 'trq_max')

    % State variables
    q = x(1:3);
    dq = x(4:6);

    % Control parameters
    K = diag([K_l K_phi]);
    D = diag([D_l D_phi]);

    switch control_mode
        case 1  % Joint Impedance
            K = diag([10 10]);
            D = diag([0.1 0.1]);

            q_cur = q(2:3);
            [q1, q2] = polar2angle(x_set(1), x_set(2), l_UL, l_LL);
            q_set = [q1; q2];
            e = q_set - q_cur;

            u = K * e - D * dq(2:3);

        case 2  % Cartesian Impedance
            K = diag([300 300]);
            D = diag([6 6]);

            J_cq = get_J_cq(q);
            dc = J_cq * dq(2:3);

            q_cur = get_foot_pos(q);
            [x_setpoint, y_setpoint] = polar2cartesian(x_set(1), x_set(2));
            q_set = [x_setpoint; y_setpoint];
            e = q_set - q_cur;

            F_c = K * e - D * dc;
            u = J_cq' * F_c;

        case 3  % Polar Impedance
            J_pq = get_J_pq(q);
            dr_p = J_pq * dq(2:3);

            x_cur = get_r_F_pol(q);
            e = x_set - x_cur;
            
            % Tutorial 5: Equation (3.3)
            F_p = K * e - D * dr_p;
            
            % Tutorial 5: Equation (3.5)
            u = J_pq' * F_p;
    end
    
    if enable_trq_limit
        u = limit_torque(u, trq_max);
    end

end

function u = limit_torque(u, trq_max)
    u(u > trq_max) = trq_max;
    u(u < -trq_max) = -trq_max;
end

function [q1,q2] = polar2angle(l, phi, l1, l2)
    q1 = phi - acos(l1/l);
    q2 = phi - q1 + acos(l2/l);
end

function [x,y] = polar2cartesian(l, phi)
    x = l * cos(phi);
    y = l * sin(phi);
end
