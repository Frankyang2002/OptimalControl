import casadi.*

% Parameters
m = 80;         % Mass
k = 2000;       % Spring stiffness
l0 = 5;         % Natural spring length
d = 0.01;       % Damping coefficient
g = 9.8;        % Gravitational acceleration
N = 100;        % Number of control intervals
T = 2;          % Fixed total time

% ---- Optimization problem --------
opti = casadi.Opti(); % Optimization problem

% Decision variables
X = opti.variable(4, N+1); % State trajectory: [x_c; y_c; dx_c; dy_c]
pos_x = X(1,:);           % Horizontal position
pos_y = X(2,:);           % Vertical position
vel_x = X(3,:);           % Horizontal velocity
vel_y = X(4,:);           % Vertical velocity
U = opti.variable(1, N);  % Control input (external force)

% ---- Objective --------
opti.minimize(sumsqr(U)); % Minimize control effort (sum of u^2)

% ---- Dynamic constraints --------
dt = T / N; % Time step size
f = @(x, u) spring_dynamics(x, u, k, l0, d, m, g); % Dynamics function

for k = 1:N
    % Runge-Kutta 4 integration
    k1 = f(X(:,k),         U(:,k));
    k2 = f(X(:,k)+dt/2*k1, U(:,k));
    k3 = f(X(:,k)+dt/2*k2, U(:,k));
    k4 = f(X(:,k)+dt*k3,   U(:,k));
    x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4); 
    opti.subject_to(X(:,k+1) == x_next); % Ensure trajectory consistency
end

% ---- Path constraints --------
opti.subject_to(pos_y >= 0);  % Vertical position must stay above ground
%opti.subject_to(-500 <= U <= 500); % Control input is limited
opti.subject_to(-2000 <= U <= 2000);

% ---- Boundary conditions --------
opti.subject_to(pos_x(1) == -3);   % Initial x position
opti.subject_to(pos_y(1) == 4);    % Start at ground level
opti.subject_to(vel_x(1) == 3);    % Initial horizontal velocity
opti.subject_to(vel_y(1) == 0);    % Initial vertical velocity
opti.subject_to(vel_x(N+1) == 3);  % Target x velocity
opti.subject_to(vel_y(N+1) == 0);  % Target y velocity
opti.subject_to(pos_x(N+1) == 3);  % Target x position
%opti.subject_to(pos_x(N+1) <= 4);
opti.subject_to(pos_y(N+1) == 4);  % Target y position
%opti.subject_to(pos_y(N+1) <= 5);

% ---- Initial guess --------
opti.set_initial(pos_y, linspace(4, 4, N+1)); % Keep y constant as an initial guess
opti.set_initial(U, 100*ones(1, N));            % Initial guess for control

% ---- Solver --------
opti.solver('ipopt'); % Set solver
sol = opti.solve();   % Solve optimization problem

% ---- Post-processing --------
t = linspace(0, T, N+1);

figure;
subplot(2,1,1); 
hold on;
plot(t, sol.value(pos_x), 'r', 'LineWidth', 1.5);
plot(t, sol.value(pos_y), 'b', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Position');
legend('Horizontal Position', 'Vertical Position');
title('Optimal Position Trajectory');
grid on;

subplot(2,1,2); 
stairs(t(1:end-1), sol.value(U), 'k', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Control Input [Force]');
legend('Control Input');
title('Optimal Control Input');
grid on;

%trajectory
figure;
plot(sol.value(pos_x), sol.value(pos_y), 'm-', 'LineWidth', 2);
hold on;
plot(sol.value(pos_x(1)), sol.value(pos_y(1)), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); 
plot(sol.value(pos_x(end)), sol.value(pos_y(end)), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); 
xlabel('Horizontal Position [x]');
ylabel('Vertical Position [y]');
title('Trajectory in XY Plane');
legend('Trajectory', 'Start Point', 'End Point');
axis equal; 
grid on;

%vector
figure;
subplot(2,1,1); 
plot(t, sol.value(vel_x), 'r-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Velocity [m/s]');
title('Velocity in X Direction');
grid on;

subplot(2,1,2);
plot(t, sol.value(vel_y), 'b-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Velocity [m/s]');
title('Velocity in Y Direction');
grid on;

% ---- Dynamics function --------
function dxdt = spring_dynamics(x, u, k, l0, d, m, g)
    x_c = x(1); y_c = x(2); dx_c = x(3); dy_c = x(4);

    % Fixed foot position
    x_f = 0; 
    y_f = 0;

    % Spring dynamics
    eps = 1e-6; % To avoid division by zero
    l = sqrt((x_c - x_f)^2 + (y_c - y_f)^2) + eps; % Spring length
    dl = (dx_c * (x_c - x_f) + dy_c * (y_c - y_f)) / l; % Spring rate of change
    F_spring = -k * (l - l0) * (l > l0); % Only act when stretched
    F_damping = -d * dl;

    % Total forces
    F_total = F_spring + F_damping + u;
    Fx = F_total * (x_c - x_f) / l; % Horizontal force
    Fy = F_total * (y_c - y_f) / l - m * g; % Vertical force

    % Dynamics
    dxdt = [dx_c; dy_c; Fx / m; Fy / m];
end