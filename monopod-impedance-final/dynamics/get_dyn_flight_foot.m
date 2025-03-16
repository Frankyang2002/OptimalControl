% function dx = get_dyn_flight_foot(x, u, m, I1, I2, g)
%     % x = [y; q1; q2; dy; dq1; dq2]
%     import casadi.*
%     dx = MX.zeros(6,1);
% 
%     % 1) positions
%     dx(1)= x(4);   % y dot
%     dx(2)= x(5);   % q1 dot
%     dx(3)= x(6);   % q2 dot
% 
%     % 2) velocity derivatives
%     dd_y  = -g;   % no contact
%     dd_q1 = u(1)/I1;
%     dd_q2 = u(2)/I2;
% 
%     dx(4)= dd_y;
%     dx(5)= dd_q1;
%     dx(6)= dd_q2;
% end