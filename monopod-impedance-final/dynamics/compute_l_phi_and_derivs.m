function [l, phi, dl, dphi] = compute_l_phi_and_derivs(x, L1, L2)
    % x = [y; q1; q2; dy; dq1; dq2]
    import casadi.*

    y    = x(1);
    q1   = x(2);
    q2   = x(3);
    dy   = x(4);
    dq1  = x(5);
    dq2  = x(6);

    % (A) foot pos in the WORLD frame
    % from your get_foot_pos:
    xFoot = 0 + L1*sin(q1) + L2*sin(q1+q2);
    yFoot = y - ( L1*cos(q1) + L2*cos(q1+q2) );

    % (B) but for polar coords, we want foot position
    % RELATIVE to the "hip" joint. 
    % In a 1D hopper, let's assume the "hip" is at [0, y].
    % Then foot rel = [ xFoot - 0; yFoot - y] = [xFoot; yFoot - y].
    % Actually yFoot - y might be "negative" if foot is below. 
    % We'll do the direct approach:
    x_rel = xFoot;        % same because hip x=0
    y_rel = yFoot - y;    % shift so that hip is at y

    % l = radial distance from hip to foot
    l = sqrt( x_rel^2 + y_rel^2 );

    % phi = angle from ? (use atan2)
    phi = atan2(y_rel, x_rel);

    % (C) We also want derivatives dl, dphi
    % dl/dq, dphi/dq => but simpler to do chain rule:
    %    dl/dt   = (dl/dx_rel)*dx_rel/dt + ...
    % You can do a small step approach for demonstration:
    % But let's do it analytically:

    % dx_rel = derivative wrt time of x_rel
    % x_rel = xFoot
    % xFoot(t) = ...
    % => dxFoot = ...
    % We'll do partial derivatives or a simpler approach:
    % For now we show the final expression:

    % "Brute force" approach: define symbolic expressions and do casadi.jacobian
    % Or do direct geometry. We'll do a quick approach:

    % dxFoot/dt = d/dt of L1*sin(q1)+ L2*sin(q1+q2)
    dxFoot = L1*cos(q1)*dq1 + L2*cos(q1+q2)*(dq1 + dq2);

    % dyFoot = derivative of  [ y - (L1*cos(q1)+ L2*cos(q1+q2)) ]
    dyFoot = dy - [ -L1*sin(q1)*dq1 - L2*sin(q1+q2)*(dq1 + dq2) ];
    % carefully note the minus sign
    % => dyFoot = dy + L1*sin(q1)*dq1 + L2*sin(q1+q2)*(dq1 + dq2)

    % So x_rel_dot = dxFoot
    %    y_rel_dot = dyFoot - dy ??? Actually we already subtracted y in yFoot
    % Actually we had yFoot = y - ...
    % => y_rel = yFoot - y = [y - (stuff)] - y = - (stuff)
    % So let's keep it simpler: We'll just do:
    x_rel_dot = dxFoot;
    y_rel_dot = dyFoot - dy; 
    % Because d/dt of (yFoot - y) = dyFoot - dy

    % Then:
    % dl = (1/l) * (x_rel*x_rel_dot + y_rel*y_rel_dot)
    % dphi = (1/(x_rel^2 + y_rel^2)) * (x_rel*y_rel_dot - y_rel*x_rel_dot)
    % (these are standard 2D geometry formula for derivative of polar coords)
    dl = (1/l)*( x_rel*x_rel_dot + y_rel*y_rel_dot );   
    dphi = (1/(l^2))*( x_rel*y_rel_dot - y_rel*x_rel_dot );
end
