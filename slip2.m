import casadi.*

% Parameters
m = 80;         % Mass
k = 2000;       % Spring stiffness
l0 = 5;         % Natural spring length
d = 0.01;       % Damping coefficient
g = 9.8;        % Gravitational acceleration
N = 100;        % Number of control intervals
T = 2;          % Total time

% ---- Optimization problem --------
opti = casadi.Opti(); % Create optimization problem

% Decision variables
X = opti.variable(4, N+1); % State trajectory: [x_c; y_c; dx_c; dy_c]
pos_x = X(1,:);           % Horizontal position
pos_y = X(2,:);           % Vertical position
vel_x = X(3,:);           % Horizontal velocity
vel_y = X(4,:);           % Vertical velocity

F = opti.variable(1, N+1); % Control input (force)
U = opti.variable(1, N);   % Rate of change of F (u = dF/dt)

% ---- Objective --------
opti.minimize(sumsqr(U)); % Minimize control effort (sum of u^2)

% ---- Dynamic constraints --------
dt = T / N; % Time step
f = @(x, F) spring_dynamics(x, F, k, l0, d, m, g); % Dynamics function

for k = 1:N
    % Runge-Kutta 4 integration for dynamics
    k1 = f(X(:,k), F(:,k));
    k2 = f(X(:,k)+dt/2*k1, F(:,k)+dt/2*U(:,k));
    k3 = f(X(:,k)+dt/2*k2, F(:,k)+dt/2*U(:,k));
    k4 = f(X(:,k)+dt*k3, F(:,k)+dt*U(:,k));
    x_next = X(:,k) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    opti.subject_to(X(:,k+1) == x_next); % Ensure trajectory consistency
    
    % Update F using its derivative u
    F_next = F(:,k) + dt * U(:,k);
    opti.subject_to(F(:,k+1) == F_next); % Consistency for F
end

% ---- Path constraints --------
opti.subject_to(pos_y >= 0);  % Vertical position must stay above ground
%opti.subject_to(-5000 <= F <= 5000); % Force limits
%opti.subject_to(0 <= U ); % Rate of change of force limits

% ---- Boundary conditions --------
opti.subject_to(pos_x(1) == -3);   % Initial x position
opti.subject_to(pos_y(1) == 4);    % Initial y position
opti.subject_to(vel_x(1) == 3);    % Initial horizontal velocity
opti.subject_to(vel_y(1) == 0);    % Initial vertical velocity
opti.subject_to(U(1) == 0);        % Initial rate of change of F (u)
opti.subject_to(U(N) == 0);        % Final rate of change of F (u)
opti.subject_to(vel_x(N+1) == 3);  % Target x velocity
opti.subject_to(vel_y(N+1) == 0);  % Target y velocity
opti.subject_to(pos_x(N+1) == 3);  % Target x position
opti.subject_to(pos_y(N+1) == 4);  % Target y position

% ---- Initial guess --------
opti.set_initial(pos_y, linspace(4, 4, N+1)); % Initial guess for y position
opti.set_initial(F, 100*ones(1, N+1));           % Initial guess for force
% ---- Solver --------
opti.solver('ipopt'); % Set solver
sol = opti.solve();   % Solve optimization problem

% ---- Post-processing --------
t = linspace(0, T, N+1);

% Plot position trajectory
figure;
subplot(2,1,1); 
stairs(t(1:end-1), sol.value(U), 'g', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Rate of Change of F [u]');
title('Rate of Change of Force (u) over Time');
legend('u');
grid on;

% Plot control input F
subplot(2,1,2); 
stairs(t, sol.value(F), 'k', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Control Input [Force]');
legend('F');
title('Optimal Control Input (Force)');
grid on;

% XY plane trajectory
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

% Velocity plots
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
function dxdt = spring_dynamics(x, F, k, l0, d, m, g)
    % Unpack states
    x_c = x(1); y_c = x(2); dx_c = x(3); dy_c = x(4);

    % Fixed foot position
    x_f = 0; 
    y_f = 0;

    % Spring dynamics
    eps = 1e-6; % To avoid division by zero
    l = sqrt((x_c - x_f)^2 + (y_c - y_f)^2) + eps; % Spring length
    dl = (dx_c * (x_c - x_f) + dy_c * (y_c - y_f)) / l; % Rate of change of spring length
    F_spring = -k * (l - l0) * (l > l0); % Active only when stretched
    F_damping = -d * dl;

    % Total forces
    F_total = F_spring + F_damping + F;
    Fx = F_total * (x_c - x_f) / l; % Horizontal force
    Fy = F_total * (y_c - y_f) / l - m * g; % Vertical force

    % Dynamics equations
    dxdt = [dx_c; dy_c; Fx / m; Fy / m];
end