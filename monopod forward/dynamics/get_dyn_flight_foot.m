function dx = get_dyn_flight_foot(x, u, m, I1, I2, g)
    % x=[x; y; q1; q2; dx; dy; dq1; dq2]
    import casadi.*
    dx = MX.zeros(8,1);

    % 1) positions
    dx(1)= x(5);   % x dot
    dx(2)= x(6);   % y dot
    dx(3)= x(7);   % q1 dot
    dx(4)= x(8);   % q2 dot

    % 2) velocity derivatives
    dd_x  = 0;
    dd_y  = -g;   % no contact
    dd_q1 = u(1)/I1;
    dd_q2 = u(2)/I2;

    dx(5)= dd_x;
    dx(6)= dd_y;
    dx(7)= dd_q1;
    dx(8)= dd_q2;
end