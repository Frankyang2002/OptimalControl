% function dx = get_dyn_stance_foot(x, u, m, I1, I2, g)
%     % stance => foot on ground => assume foot y=0, so reaction force
%     import casadi.*
%     dx = MX.zeros(6,1);
% 
%     dx(1)= x(4);
%     dx(2)= x(5);
%     dx(3)= x(6);
% 
%     k_leg=2000;
%     F_leg = -k_leg * ( x(1) - 0 );
% 
%     dd_y = F_leg/m - g;
%     dd_q1= u(1)/I1;
%     dd_q2= u(2)/I2;
% 
%     dx(4)= dd_y;
%     dx(5)= dd_q1;
%     dx(6)= dd_q2;
% end