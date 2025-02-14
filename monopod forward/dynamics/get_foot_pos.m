function p_foot = get_foot_pos(x, L1, L2)
    % x=[x; y; q1; q2; dx; dy; dq1; dq2]
    xpos   = x(1);
    ypos   = x(2);
    q1  = x(3);
    q2  = x(4);
    xFoot = xpos + L1*cos(q1) + L2*cos(q1+q2-pi);  
    yFoot = ypos - L1*sin(q1) - L2*sin(q1+q2-pi);
    p_foot= [xFoot; yFoot];
end