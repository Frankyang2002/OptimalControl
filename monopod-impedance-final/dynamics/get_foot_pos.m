function p_foot = get_foot_pos(x, L1, L2)
    % x=[y; q1; q2; dy; dq1; dq2]
    y   = x(1);
    q1  = x(2);
    q2  = x(3);
    xFoot = 0 + L1*sin(q1) + L2*sin(q1+q2);  
    yFoot = y - ( L1*cos(q1) + L2*cos(q1+q2) );
    p_foot= [xFoot; yFoot];
end