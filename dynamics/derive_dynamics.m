%% Derive monopod dynamics
% Init
clear; clc;

syms y q1 q2 real
syms dy dq1 dq2 real
syms u1 u2 real
syms H_d
syms t

%% State & control vectors
q = [y; q1; q2];
dq = [dy; dq1; dq2];
u = [0; u1; u2];

%% Load monopod parameters
load('monopod_parameters', 'm_B', 'm_UL', 'm_LL', 'l_B', 'l_UL', 'l_LL', 'I_B', 'I_UL', 'I_LL', 'K_l', 'q_rest')
g = 9.81;

%% Homogeneous transformations
% World to base frame
W_X_B = [1, 0, 0, 0;
         0, 1, 0, y;
         0, 0, 1, 0;
         0, 0, 0, 1];

% Base to upper leg frame
B_X_UL = [cos(q1), -sin(q1), 0, l_UL/2*cos(q1);
          sin(q1),  cos(q1), 0, l_UL/2*sin(q1);
                0,        0, 1,              0;
                0,        0, 0,              1];

% Upper leg to lower leg frame
UL_X_LL = [cos(q2), -sin(q2), 0, l_UL/2+l_LL/2*cos(q2);
           sin(q2),  cos(q2), 0,        l_LL/2*sin(q2);
                 0,        0, 1,                     0;
                 0,        0, 0,                     1];

% Lower leg to foot frame
LL_X_F = [1, 0, 0, l_LL/2;
          0, 1, 0,      0;
          0, 0, 1,      0;
          0, 0, 0,      1];

% Compute frames
W_X_UL = W_X_B * B_X_UL;
W_X_LL = W_X_UL * UL_X_LL;
W_X_F = W_X_LL * LL_X_F;
B_X_LL = B_X_UL * UL_X_LL;
B_X_F = B_X_LL * LL_X_F;

X_CoM = {W_X_B, W_X_UL, W_X_LL};

%% Polar coordinates
B_x_F = B_X_F(1, 4);
B_y_F = B_X_F(2, 4);

l_pol = sqrt(B_x_F^2 + B_y_F^2);  % Tutorial 5: Equation (3.1)
phi_pol = atan2(B_y_F, B_x_F);  % Tutorial 5: Equation (3.2)

r_F_pol = simplify([l_pol; phi_pol]);

J_pq = simplify(jacobian(r_F_pol, [q1, q2]));  % Polar jacobian

%% Translational jacobian
for i = 1:length(X_CoM)
    r_CoM{i} = simplify(X_CoM{i}(1:3, 4));  % Position vector
    J_v{i} = jacobian(r_CoM{i}, q);
    J_w{i} = sym(zeros(3,3));
end

% Manually derived from rotational joints
J_w{2}(3,3) = 1;
J_w{3}(3,3) = 1;

%% Mass matrix
m = {m_B, m_UL, m_LL};  % Mass
I = {I_B, I_UL, I_LL};  % Inertia
M = sym(zeros(3, 3));

for i = 1:length(m)
    M = M + m{i} * (transpose(J_v{i}) * J_v{i}) + (transpose(J_w{i}) * I{i} * J_w{i});  % Tutorial 1: Equation 4.15
end

M = simplify(M);
M_inv = inv(M);

%% Coriolis Matrix
C = sym(zeros(3, 3));

for i=1:size(C,1)
    for j=1:size(C,2)
        for k=1:length(q)
           C(i,j) = C(i,j) + 0.5 * (diff(M(i, j), q(k)) + diff(M(i, k), q(j)) - diff(M(j, k), q(i))) * dq(i);  % Tutorial 1: Equation 4.22
        end
    end
end

C = simplify(C);

%% Gravitational terms
grav = [0; g; 0];
E_pot_g = sym(0.);

for i = 1:length(X_CoM)
    E_pot_g = E_pot_g + (m{i} * r_CoM{i});
end
E_pot_g = transpose(grav) * E_pot_g;

G = simplify(transpose(jacobian(E_pot_g, q)));  % Tutorial 1 Equation 4.24

%% Rest length for elastic potential energy
q1_rest = q_rest(2);
q2_rest = q_rest(3);
l_rest = subs(l_pol, [q1 q2], [q1_rest q2_rest]);

%% Total system energy
E_pot_e = 0.5 * K_l * (l_rest - l_pol)^2;
E_kin = 0.5 * dq' * M * dq;

E_sys = E_pot_g + E_kin + E_pot_e;

%% Desired system energy at high
B_y_UL = B_X_UL(2, 4);
B_y_LL = B_X_LL(2, 4);

F_y_UL = l_rest + B_y_UL;  % l_rest = -B_y_F (y-axis in B frame points upwards)
F_y_LL = l_rest + B_y_LL;

y_B_d = H_d + l_rest;
y_UL_d = H_d + F_y_UL;
y_LL_d = H_d + F_y_LL;

y_d = [y_B_d y_UL_d y_LL_d];

E_sys_d = [m_B m_UL m_LL] * g * y_d';

%% Constraint forces
r_F = simplify(W_X_F(1:3, 4));  % Position vector of foot
J_r_F = jacobian(r_F, q);  % Jacobian of foot position
J_c = simplify(J_r_F(1:2, :));  % Constraint Jacobian

for i = 1:2
    H{i} = jacobian(J_c(i,:), q);  % Hessian tensor
    dJ_c{i} = H{i} * dq;
end

J_c_T = transpose(J_c);
J_c_dot = simplify([transpose(dJ_c{1}); transpose(dJ_c{2})]);

n = C * dq + G;

L = simplify(inv(J_c * M_inv * J_c_T));
lambda = L * (J_c * M_inv * (n - u) - J_c_dot * dq);  % Sobotka: Equation 2.15
lambda = simplify(lambda);

%% Impact law
dq_plus = dq - M_inv * J_c_T * L * J_c * dq;  % Sobotka: Equation 2.19
dq_plus = simplify(dq_plus);

%% Flight phase
ddq_flight = M_inv * (u - n);  % Tutorial 5: Equation 1.1 solved for ddq

%% Stance phase
ddq_stance = M_inv * (u + J_c_T * lambda - n);  % Tutorial 5: Equation 1.2 solved for ddq

%% ODEs
ode_flight = [dq; ddq_flight];
ode_stance = [dq; ddq_stance];

%% Export functions
if isfolder('dynamics')
    path_dest = fullfile(pwd, 'dynamics');
else
    path_dest = fullfile(fileparts(pwd), 'dynamics');
end

% Foot y position wrt world frame
W_y_F = W_X_F(2, 4);
matlabFunction(W_y_F, 'File', fullfile(path_dest, 'get_c'), 'Vars', {q});

% Foot position wrt body frame
foot_pos = [B_X_F(1:2, 4)];
matlabFunction(foot_pos, 'File', fullfile(path_dest, 'get_foot_pos'), 'Vars', {q});

% Frame orientation and position wrt world frame
matlabFunction(W_X_B, 'File', fullfile(path_dest, 'get_W_X_B'), 'Vars', {q});
matlabFunction(W_X_UL, 'File', fullfile(path_dest, 'get_W_X_UL'), 'Vars', {q});
matlabFunction(W_X_LL, 'File', fullfile(path_dest, 'get_W_X_LL'), 'Vars', {q});

% Foot position in polar copordinates wrt world frame
matlabFunction(r_F_pol, 'File', fullfile(path_dest, 'get_r_F_pol'), 'Vars', {q});

% Jacobian polar to general coordinates
matlabFunction(J_pq, 'File', fullfile(path_dest, 'get_J_pq'), 'Vars', {q});

% Jacobian cartesian to general coordinates
J_cq = simplify(jacobian([B_x_F,B_y_F],[q1,q2]));
matlabFunction(J_cq, 'File', fullfile(path_dest, 'get_J_cq'), 'Vars', {q});

% Jacobians
matlabFunction(J_v{1}, 'File', fullfile(path_dest, 'get_J_v_B'));
matlabFunction(J_v{2}, 'File', fullfile(path_dest, 'get_J_v_UL'), 'Vars', {q});
matlabFunction(J_v{3}, 'File', fullfile(path_dest, 'get_J_v_LL'), 'Vars', {q});

% Rest length for elastic potential energy
matlabFunction(l_rest, 'File', fullfile(path_dest, 'get_l_rest'));

% Mass matrix
matlabFunction(M, 'File', fullfile(path_dest, 'get_M'), 'Vars', {q});

% Bias term
matlabFunction(n, 'File', fullfile(path_dest, 'get_n'), 'Vars', {q, dq});

% Contact Jacobian
matlabFunction(J_c, 'File', fullfile(path_dest, 'get_J_c'), 'Vars', {q});


% System energies
matlabFunction(E_sys, 'File', fullfile(path_dest, 'get_E_sys'), 'Vars', {[q; dq]});
matlabFunction(E_sys_d, 'File', fullfile(path_dest, 'get_E_sys_d'), 'Vars', {q; H_d});

% Contact forces function
matlabFunction(lambda, 'File', fullfile(path_dest, 'get_lambda'), 'Vars', {[q; dq]; [u1; u2]});

% Flight dynamics function
matlabFunction(ode_flight, 'File', fullfile(path_dest, 'get_dyn_flight'), 'Vars', {t; [q; dq]; [u1; u2]});

% Stance dynamics function
matlabFunction(ode_stance, 'File', fullfile(path_dest, 'get_dyn_stance'), 'Vars', {t; [q; dq]; [u1; u2]});

% Impact jump map function
matlabFunction(dq_plus, 'File', fullfile(path_dest, 'get_impact'), 'Vars', {[q; dq]});
