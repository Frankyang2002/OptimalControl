function dx = get_dyn(t, x)

q = x(1:3);
dq = x(4:6);

lambda = get_contact_force(q,dq);
J_c = get_J_c(q);
J_c = J_c(2,:);
M = get_M(q);
h = get_n(q,dq);

% Dynamics
%M * ddq + h + J_c' * lambda = 0;

dx = [dq;
      M\ (-h - J_c' * lambda)];


end