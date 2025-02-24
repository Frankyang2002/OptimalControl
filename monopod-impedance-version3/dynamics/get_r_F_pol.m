function r_F_pol = get_r_F_pol(x)
    % Compute foot position in polar coordinates
    q1 = x(2);
    q2 = x(3);
    
    % Link lengths
    L1 = 0.5;
    L2 = 0.5;
    
    % Compute foot position
    xFoot = L1 * sin(q1) + L2 * sin(q1 + q2);
    yFoot = - (L1 * cos(q1) + L2 * cos(q1 + q2));
    
    % Convert to polar coordinates
    l = sqrt(xFoot^2 + yFoot^2);
    phi = atan2(yFoot, xFoot);
    
    r_F_pol = [l; phi];
end
