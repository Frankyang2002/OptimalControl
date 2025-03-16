function Jp = get_polar_jacobian(x, L1, L2)
    import casadi.*

    % x=[y; q1; q2; dy; dq1; dq2]
    q1 = x(2);
    q2 = x(3);

    % As in compute_l_phi_and_derivs:
    % But let's define them here again for clarity
    x_rel = L1*sin(q1) + L2*sin(q1 + q2);
    y_rel = - (L1*cos(q1) + L2*cos(q1 + q2));  
    % so y_rel = -(...)= -L1*cos(q1) - L2*cos(q1+q2)
    
    % l, phi
    l    = sqrt(x_rel^2 + y_rel^2);
    phi  = atan2(y_rel, x_rel);

    %-------------------------------------------
    % Compute partial derivatives by hand
    %-------------------------------------------
    % We'll define short-hands:
    %  dx_rel/dq1, dx_rel/dq2, dy_rel/dq1, dy_rel/dq2:

    dx_rel_dq1 = L1*cos(q1) + L2*cos(q1 + q2);
    dx_rel_dq2 = L2*cos(q1 + q2);

    % y_rel = - [L1*cos(q1) + L2*cos(q1+q2)]
    % => dy_rel/dq1 = - [ -L1*sin(q1) - L2*sin(q1+q2) ] = L1*sin(q1) + L2*sin(q1+q2)
    dy_rel_dq1 = L1*sin(q1) + L2*sin(q1 + q2);
    dy_rel_dq2 = L2*sin(q1 + q2);

    %-------------------------------------------
    %  dl/dq1 and dl/dq2
    %-------------------------------------------
    %  l = sqrt(x_rel^2 + y_rel^2)
    %  dl/dq1 = (1/l) * [ x_rel*(dx_rel/dq1) + y_rel*(dy_rel/dq1) ]
    
    dl_dq1 = (1/l)*( x_rel*dx_rel_dq1 + y_rel*dy_rel_dq1 );
    dl_dq2 = (1/l)*( x_rel*dx_rel_dq2 + y_rel*dy_rel_dq2 );

    %-------------------------------------------
    %  dphi/dq1, dphi/dq2
    %-------------------------------------------
    %  phi = atan2(y_rel, x_rel)
    %  dphi/dq1 = 1/(x_rel^2 + y_rel^2) * [ x_rel*(dy_rel/dq1) - y_rel*(dx_rel/dq1) ]
    denom = (x_rel^2 + y_rel^2);

    dphi_dq1 = (1/denom)*( x_rel*dy_rel_dq1 - y_rel*dx_rel_dq1 );
    dphi_dq2 = (1/denom)*( x_rel*dy_rel_dq2 - y_rel*dx_rel_dq2 );

    %-------------------------------------------
    % Construct the Jacobian
    %-------------------------------------------
    %    [ dl/dq1    dl/dq2   ]
    % J = [                     ]
    %    [ dphi/dq1  dphi/dq2 ]
    %
    Jp = MX(2,2);
    Jp(1,1) = dl_dq1;
    Jp(1,2) = dl_dq2;
    Jp(2,1) = dphi_dq1;
    Jp(2,2) = dphi_dq2;
end
