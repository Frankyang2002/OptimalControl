function dx = get_dyn_stance_foot(x, u, m, I1, I2, g, L1, L2, LRest, k_leg)
    % x=[x; y; q1; q2; dx; dy; dq1; dq2]
    % stance => foot on ground => assume foot y=0
    import casadi.*
    dx = MX.zeros(8,1);
    
    x_hip = x(1);
    y_hip = x(2);
    foot_pos = get_foot_pos(x,L1,L2);
    x_foot = foot_pos(1);
    y_foot = foot_pos(2);
    
    % Compute relative position of hip to foot
    dx_hip_foot = x_hip - x_foot;  
    dy_hip_foot = y_hip - y_foot;  

    % Simplify as a spring model
    L = sqrt((dx_hip_foot)^2 + (dy_hip_foot)^2);
    F_leg = k_leg * (LRest - L);
    F_angle = atan2(dy_hip_foot, dx_hip_foot);

    % Pushes up the leg
    Fx = F_leg * cos(F_angle);
    Fy = F_leg * sin(F_angle);  

    dx(1)= x(5);
    dx(2)= x(6); 
    dx(3)= x(7);
    dx(4)= x(8);

    dd_x = (Fx / m);
    dd_y = (Fy / m) - g;
    dd_q1= u(1) / I1;
    dd_q2= u(2) / I2;

    dx(5)= dd_x;
    dx(6)= dd_y;
    dx(7)= dd_q1;
    dx(8)= dd_q2;
    
end