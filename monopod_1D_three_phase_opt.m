%function monopod_1D_three_phase_opt()
    import casadi.*;

    % Create Opti instance
    opti = Opti();

    % --- 1) Decision variables for phase durations ---
    T_f1 = opti.variable();  % Time for flight down
    T_s  = opti.variable();  % Time for stance
    T_f2 = opti.variable();  % Time for flight up

    % Enforce minimal positive time
    opti.subject_to(T_f1 >= 0);
    opti.subject_to(T_s  >= 0);
    opti.subject_to(T_f2 >= 0);
    opti.subject_to(T_f1 + T_s + T_f2 <= 15)

    % Fixed number of intervals for each phase
    N1 = 20;
    N2 = 20;
    N3 = 20;

    % Symbolic time steps
    dt_f1 = T_f1 / N1;
    dt_s  = T_s  / N2;
    dt_f2 = T_f2 / N3;

    % --- 2) Physical parameters ---
    g  = 9.81;
    m  = 10.0;   % mass
    I1 = 0.02;   % inertia joint1
    I2 = 0.01;   % inertia joint2

    % --- 3) States and inputs ---
    % x = [y; q1; q2; dy; dq1; dq2]
    Xf1 = opti.variable(6, N1+1);  
    Uf1 = opti.variable(2, N1);   

    Xs  = opti.variable(6, N2+1);  
    Us  = opti.variable(2, N2);

    Xf2 = opti.variable(6, N3+1);  
    Uf2 = opti.variable(2, N3);

    % --- 4) Objective: sum of squared torques + soft periodicity ---
    alpha = 100000;  % large weight for soft periodic condition
    % e = (final state) - (initial state)
    e = Xf2(:,N3+1) - Xf1(:,1);

    obj = sumsqr(Uf1) + sumsqr(Us) + sumsqr(Uf2) ...
          + alpha * sumsqr(e);
    opti.minimize(obj);

    % --- 5) Dynamics constraints (RK4) ---
    % Flight dynamics
    f_flight = @(x,u) get_dyn_flight_yOnly(x,u,m,g,I1,I2);

    for k = 1:N1
        xk  = Xf1(:, k);
        uk  = Uf1(:, k);
        xk1 = rk4_step(f_flight, xk, uk, dt_f1);
        opti.subject_to(Xf1(:, k+1) == xk1);
    end

    for k = 1:N3
        xk  = Xf2(:, k);
        uk  = Uf2(:, k);
        xk1 = rk4_step(f_flight, xk, uk, dt_f2);
        opti.subject_to(Xf2(:, k+1) == xk1);
    end

    % Stance dynamics
    f_stance = @(x,u) get_dyn_stance_yOnly(x,u,m,g,I1,I2);

    for k = 1:N2
        xk  = Xs(:, k);
        uk  = Us(:, k);
        xk1 = rk4_step(f_stance, xk, uk, dt_s);
        opti.subject_to(Xs(:, k+1) == xk1);
    end

    % --- 6) Initial conditions ---
    % At t=0, y=1, dy=0, q1=-135deg, q2=90deg, dq=0
    opti.subject_to(Xf1(1,1) == 1);
    opti.subject_to(Xf1(4,1) == 0);
    opti.subject_to(Xf1(2,1) == deg2rad(-135));
    opti.subject_to(Xf1(3,1) == deg2rad(90));
    opti.subject_to(Xf1(5,1) == 0);
    opti.subject_to(Xf1(6,1) == 0);

    % Phase1 -> Phase2 continuity
    opti.subject_to(Xf1(1,N1+1) == 0); % get_C(xf(end)) = 0
    opti.subject_to( Xs(:,1) == Xf1(:,N1+1) ); % Xs == getImpact(xf)
    % Phase2 -> Phase3 continuity
    %opti.subject_to(Xs(1,N2+1) == 0); 
    %get_lambda(xs, us) = 0
    opti.subject_to( Xf2(:,1) == Xs(:,N2+1) );

    % --- 7) Path constraints ---
    % y >= 0 in all phases
    opti.subject_to(Xf1(1,:) >= 0);
    opti.subject_to(Xs(1,:)  >= 0);
    opti.subject_to(Xf2(1,:) >= 0);

    % Torque limits
    maxTorque = 20000;
    opti.subject_to( -maxTorque <= Uf1(:) <= maxTorque );
    opti.subject_to( -maxTorque <= Us(:)  <= maxTorque );
    opti.subject_to( -maxTorque <= Uf2(:) <= maxTorque );

    % --- 8) Initial guesses ---
    opti.set_initial(Xf1, 0);
    opti.set_initial(Xs,  0);
    opti.set_initial(Xf2, 0);
    opti.set_initial(Uf1, 0);
    opti.set_initial(Us,  0);
    opti.set_initial(Uf2, 0);

    % Provide initial guesses for times
    opti.set_initial(T_f1, 1);
    opti.set_initial(T_s,  1);
    opti.set_initial(T_f2, 1);

    % --- 9) Solve ---
    opti.solver('ipopt');
    sol = opti.solve();
%%
    % Retrieve final solution
    T_f1_sol =opti.debug.value(T_f1);
    T_s_sol  =opti.debug.value(T_s);
    T_f2_sol =opti.debug.value(T_f2);

    disp(['Optimal T_f1=', num2str(T_f1_sol), ...
          ', T_s=', num2str(T_s_sol), ...
          ', T_f2=', num2str(T_f2_sol)]);

    Xf1_sol =opti.debug.value(Xf1);
    Xs_sol  =opti.debug.value(Xs);
    Xf2_sol =opti.debug.value(Xf2);
    Uf1_sol =opti.debug.value(Uf1);
    Us_sol  =opti.debug.value(Us);
    Uf2_sol =opti.debug.value(Uf2);

  %% --- 10) Visualization ---

% Build time vectors for the states (which have N+1 nodes each)
tf1_vec = linspace(0, T_f1_sol, N1+1);
ts_vec  = linspace(T_f1_sol, T_f1_sol + T_s_sol, N2+1);
tf2_vec = linspace(T_f1_sol + T_s_sol, T_f1_sol + T_s_sol + T_f2_sol, N3+1);

figure('Name','Monopod Y Trajectory','Color','white');
hold on; grid on;
plot(tf1_vec, Xf1_sol(1,:), 'r-o','DisplayName','Flight Down');
plot(ts_vec,  Xs_sol(1,:),  'g-o','DisplayName','Stance');
plot(tf2_vec, Xf2_sol(1,:), 'b-o','DisplayName','Flight Up');
xlabel('Time (s)');
ylabel('y (m)');
legend();
title('Vertical Position (y)');

%% --- 10.1) Plot torque vs. time for each phase ---

% Each phase has N control intervals, so we create time vectors of length N
t_f1 = linspace(0, T_f1_sol, N1);                     % flight down
t_s  = linspace(T_f1_sol, T_f1_sol + T_s_sol, N2);    % stance
t_f2 = linspace(T_f1_sol + T_s_sol, ...
                T_f1_sol + T_s_sol + T_f2_sol, N3);   % flight up

figure('Name','Monopod Control Torques','Color','white');

% (a) Plot u1 in a subplot
subplot(2,1,1);
hold on; grid on;
plot(t_f1, Uf1_sol(1,:), 'r', 'LineWidth',1.5);
plot(t_s,  Us_sol(1,:),  'g', 'LineWidth',1.5);
plot(t_f2, Uf2_sol(1,:), 'b', 'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Torque u_1 (Nm)');
legend('Flight Down','Stance','Flight Up','Location','best');
title('Torque u_1 vs. Time');

% (b) Plot u2 in another subplot
subplot(2,1,2);
hold on; grid on;
plot(t_f1, Uf1_sol(2,:), 'r', 'LineWidth',1.5);
plot(t_s,  Us_sol(2,:),  'g', 'LineWidth',1.5);
plot(t_f2, Uf2_sol(2,:), 'b', 'LineWidth',1.5);
xlabel('Time (s)');
ylabel('Torque u_2 (Nm)');
legend('Flight Down','Stance','Flight Up','Location','best');
title('Torque u_2 vs. Time');


%% ================== RK4 one-step ==================
function x_next = rk4_step(f, x, u, dt)
    k1 = f(x, u);
    k2 = f(x + dt/2*k1, u);
    k3 = f(x + dt/2*k2, u);
    k4 = f(x + dt*k3,   u);
    x_next = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

%% ================== Flight dynamics ==================
function dx = get_dyn_flight_yOnly(x, u, m, g, I1, I2)
    % x = [y; q1; q2; dy; dq1; dq2]
    dx = [
        x(4);   % dy
        x(5);   % dq1
        x(6);   % dq2
       -g;      % ddy
        u(1)/I1;% ddq1
        u(2)/I2 % ddq2
    ];
end

%% ================== Stance dynamics ==================
function dx = get_dyn_stance_yOnly(x, u, m, g, I1, I2)
    % Simple spring force anchored at y=0, plus gravity and torque
    k_leg = 1000;
    y  = x(1);
    dy = x(4);
    F_leg = -k_leg*(y - 0);

    dx = [
        x(4);          % dy
        x(5);          % dq1
        x(6);          % dq2
        (F_leg/m)- g;  % ddy
        u(1)/I1;       % ddq1
        u(2)/I2        % ddq2
    ];
end