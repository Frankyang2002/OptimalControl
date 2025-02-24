function p_foot = get_foot_pos(x, L1, L2)
    % Compute foot position
    y = x(1);
    q1 = x(2);
    q2 = x(3);

    xFoot = L1*sin(q1) + L2*sin(q1+q2);  
    yFoot = y - ( L1*cos(q1) + L2*cos(q1+q2) );
    
    % deleted to try yFoot = max(min(yFoot, 10), -10);
    
    p_foot= [xFoot; yFoot];

    % Debugging: Display foot position
    disp('Checking Foot Position:');
    disp(p_foot);
end
