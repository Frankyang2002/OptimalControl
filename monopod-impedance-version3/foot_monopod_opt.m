function foot_monopod_opt()
    import casadi.*;
    
    opti = casadi.Opti();
    
    % !minjie! Define optimization variables for impedance control:
    l0_f1 = opti.variable();  % Desired leg length in flight phase 1
    phi0_f1 = opti.variable(); % Desired hip angle in flight phase 1
    l0_s  = opti.variable();  % Desired leg length in stance phase
    phi0_s = opti.variable(); % Desired hip angle in stance phase
    l0_f2 = opti.variable();  % Desired leg length in flight phase 2
    phi0_f2 = opti.variable(); % Desired hip angle in flight phase 2
    
    % minjie Constraints on l0
    opti.subject_to(0.8 <= l0_f1);
    opti.subject_to(l0_f1 <= 1.2);

    opti.subject_to(0.5 <= l0_s);
    opti.subject_to(l0_s <= 1.5);

    opti.subject_to(0.8 <= l0_f2);
    opti.subject_to(l0_f2 <= 1.2);

    % Constraints on phi0
    opti.subject_to(deg2rad(-60) <= phi0_f1);
    opti.subject_to(phi0_f1 <= deg2rad(60));

    opti.subject_to(deg2rad(-90) <= phi0_s);
    opti.subject_to(phi0_s <= deg2rad(90));

    opti.subject_to(deg2rad(-90) <= phi0_f2);
    opti.subject_to(phi0_f2 <= deg2rad(45));





    % 3-phase times: T_f1, T_s, T_f2
    T_f1 = opti.variable();  % flight down
    T_s  = opti.variable();  % stance
    T_f2 = opti.variable();  % flight up
    opti.subject_to(T_f1 >= 0.2);
    opti.subject_to(T_s  >= 0.2); %minjie avoid negative value
    opti.subject_to(T_f2 >= 0.2);
    opti.subject_to(T_f1 + T_s + T_f2 <= 15);
    opti.subject_to(T_f1 + T_s + T_f2 >= 1.0); % ensure time not too short

    
    % minjie
    opti.subject_to(T_f1 <= 5);
    opti.subject_to(T_s <= 5);
    opti.subject_to(T_f2 <= 5);


    % number of intervals each phase
    N1 = 20; N2 = 20; N3 = 20;
    dt_f1 = T_f1 / N1;
    dt_s  = T_s  / N2;
    dt_f2 = T_f2 / N3;

    %Model parameters 
    m = 10;   % body mass
    I1    = 0.02; % inertia of joint1
    I2    = 0.01; % inertia of joint2
    g     = 9.81; 
    L1    = 0.5;  % link1 length
    L2    = 0.5;  % link2 length
    
    % Impedance control parameters (modify these values as needed)
    K_l = 150;     % Stiffness for leg length control
    D_l = 3;       % Damping for leg length control
    K_phi = 80;    % Stiffness for hip angle control
    D_phi = 1;      % Damping for hip angle control


    %Decision variables: Xf1, Xs, Xf2; Uf1, Us, Uf2
    % State vector x = [y; q1; q2; dy; dq1; dq2], dimension=6
    Xf1 = opti.variable(6, N1+1);

    
    % minjie limit initial y position
    opti.subject_to( 0.5 <= Xf1(1,1));
    opti.subject_to( Xf1(1,1) <= 1.5);

    Xs  = opti.variable(6, N2+1);


    Xf2 = opti.variable(6, N3+1);


    %Objective function: sum of torque^2 + soft periodic
    slack = opti.variable(6,1); 
    opti.subject_to(slack >= 0);
    alpha = 50;
    e = Xf2(:, N3+1) - Xf1(:,1);  % final state - initial state

    obj= 1e-3 * sumsqr(Xf1) + 1e-3 * sumsqr(Xs) + 1e-3 * sumsqr(Xf2) ...
              + 500 * sumsqr(l0_f1 - 1.0) + 500 * sumsqr(l0_s - 1.0) + 500 * sumsqr(l0_f2 - 1.0) ...
              + 200 * sumsqr(phi0_f1) + 200 * sumsqr(phi0_s) + 200 * sumsqr(phi0_f2) ...
              + 1e3 * sumsqr(slack);
    opti.minimize(obj);

    %Build RK4 constraints for each phase
    % flight down

    for k=1:5:N1
        xk = Xf1(:,k);
        xk1 = rk4_step(@(xx) get_dyn_flight_foot(xx, l0_f1, phi0_f1, m, I1, I2, g, K_l, D_l, K_phi, D_phi), xk, dt_f1);
        opti.subject_to(Xf1(:,k+1) == xk1);
    end

    % stance

    for k=1:N2
        xk = Xs(:,k);
        xk1 = rk4_step(@(xx) get_dyn_stance_foot(xx, l0_s, phi0_s, m, I1, I2, g, K_l, D_l, K_phi, D_phi), xk, dt_s);

        opti.subject_to(Xs(:,k+1) == xk1);
    end

    % flight up

    for k=1:N3
        xk = Xf2(:,k);
        xk1 = rk4_step(@(xx) get_dyn_flight_foot(xx, l0_f2, phi0_f2, m, I1, I2, g, K_l, D_l, K_phi, D_phi), xk, dt_f2);

        opti.subject_to(Xf2(:,k+1) == xk1);
    end

    % Boundary conditions 
    % initial condition: 
    opti.subject_to(0.5 <= Xf1(1,1)); % y 
    opti.subject_to(Xf1(1,1) <= 1.5); % y 
    opti.subject_to(-0.5 <= Xf1(4,1)); % dy 
    opti.subject_to(Xf1(4,1) <= 0.5); % dy 
    opti.subject_to(deg2rad(-60) <= Xf1(2,1)); % q1 
    opti.subject_to(Xf1(2,1) <= deg2rad(-30)); % q1 
    opti.subject_to(deg2rad(60) <= Xf1(3,1)); % q2 
    opti.subject_to(Xf1(3,1) <= deg2rad(120)); % q2 
    opti.subject_to(-0.5 <= Xf1(5,1)); % dq1
    opti.subject_to(Xf1(5,1) <= 0.5); % dq1
    opti.subject_to(-0.5 <= Xf1(6,1)); % dq2
    opti.subject_to(Xf1(6,1) <= 0.5); % dq2

    %opti.subject_to(Xf1(5,1) >= 0);
    %opti.subject_to(Xf1(6,1) >= 0);

    % Phase1 -> Phase2 continuity at foot collision (foot y=0)
    pFoot_f1_end = get_foot_pos(Xf1(:,N1+1), L1, L2);
% minjie deleted    minjie deleted opti.subject_to( pFoot_f1_end(2) == 0 );  % touch ground
    opti.subject_to(-0.05 <= pFoot_f1_end(2));
    opti.subject_to(pFoot_f1_end(2) <= 0.05);


    % continuity in states
%  minjie deleted   opti.subject_to( Xs(:,1) == Xf1(:,N1+1));

    % Phase2 -> Phase3 continuity => foot y>0 => leaving ground
%  minjie deleted   opti.subject_to( Xf2(:,1) == Xs(:,N2+1));

% minjie deleted    opti.minimize(1000 * sumsqr(Xs(:,1) - Xf1(:,N1+1)) + 1000 * sumsqr(Xf2(:,1) - Xs(:,N2+1)));
    opti.subject_to(-0.2 <= Xs(:,1) - Xf1(:,N1+1));
    opti.subject_to(Xs(:,1) - Xf1(:,N1+1) <= 0.2);
    opti.subject_to(-0.2 <= Xf2(:,1) - Xs(:,N2+1));
    opti.subject_to(Xf2(:,1) - Xs(:,N2+1) <= 0.2);

    % minjie introduce slack variable to deal with continuity


    opti.subject_to(Xs(:,1) - Xf1(:,N1+1) <= slack); 
    opti.subject_to(Xf2(:,1) - Xs(:,N2+1) <= slack); 

%     obj = obj + 1e5 * sumsqr(slack);
%     opti.minimize(obj);


    
    %path constraints: y >= 0 => body doesn't go below ground?
    for k=1:(N1+1)
       footPos_k = get_foot_pos(Xf1(:,k), L1, L2);
       opti.subject_to( footPos_k(2) >= -0.01 );  
       %opti.subject_to( Xf1(1,k) >= 0 );
    end
    for k=1:(N3+1)
       footPos_k = get_foot_pos(Xf2(:,k), L1, L2);
       opti.subject_to( footPos_k(2) >= 0 );
       %opti.subject_to( Xf2(1,k) >= 0 );
    end
    % In stance, footPos=0 => keep foot locked? 
    for k=1:(N2+1)
       footPos_k = get_foot_pos(Xs(:,k), L1, L2);
% minjie deleted       opti.subject_to( footPos_k(2) == 0 );
        opti.subject_to(-0.1 <= footPos_k(2));
        opti.subject_to(footPos_k(2) <= 0.1);

    end



    %initial guess
%minjie deleted      opti.set_initial(Xf1, 0);
%minjie deleted     opti.set_initial(Xs,  0);
%minjie deleted     opti.set_initial(Xf2, 0);
    % Linear interpolation
    for k = 1:N1+1
        opti.set_initial(Xf1(:,k), [1 - 0.1*k/N1; deg2rad(-45 + 15*k/N1); deg2rad(90 - 10*k/N1); 0; 0; 0]);
    end

    for k = 1:N2+1
        opti.set_initial(Xs(:,k), [0.9 + 0.1*k/N2; deg2rad(-30 + 15*k/N2); deg2rad(100 - 5*k/N2); 0; 0; 0]);
    end

    for k = 1:N3+1
        opti.set_initial(Xf2(:,k), [1 - 0.1*k/N3; deg2rad(-60 + 20*k/N3); deg2rad(120 - 10*k/N3); 0; 0; 0]);
    end
    
    %minjie check if the initial values are reasonable
    disp('Checking initial values:');
    disp(['Xf1(1,1): ', num2str(1 - 0.1*1/20)]);
    disp(['Xf1(2,1): ', num2str(deg2rad(-45 + 15*1/20))]);
    disp(['Xf1(3,1): ', num2str(deg2rad(90 - 10*1/20))]);
    disp(['Xf1(4,1): ', num2str(0)]);
    disp(['Xf1(5,1): ', num2str(0)]);
    disp(['Xf1(6,1): ', num2str(0)]);


    opti.set_initial(T_f1, 1.0);
    opti.set_initial(T_s,  1.0);
    opti.set_initial(T_f2, 1.0);
    
    % minjie: Add reasonable initial values
    opti.set_initial(l0_f1, 1.0);
    opti.set_initial(phi0_f1, deg2rad(0));
    opti.set_initial(l0_s, 1.0);
    opti.set_initial(phi0_s, deg2rad(0));
    opti.set_initial(l0_f2, 1.0);
    opti.set_initial(phi0_f2, deg2rad(0));



   

    
    %Solve
%minjie     opti.solver('ipopt');
    opts = struct;
    opts.ipopt.max_iter = 5000; 
    opts.ipopt.tol = 1e-4; 
    opts.ipopt.constr_viol_tol = 1e-3; 
    opts.ipopt.dual_inf_tol = 1e-3; 
    opts.ipopt.acceptable_tol = 1e-3; 
    opts.ipopt.acceptable_constr_viol_tol = 1e-2;

    opti.solver('ipopt', opts);

%    minjie dele sol = opti.solve();
    % Solve with debugging
    try
        sol = opti.solve();
        disp('l0_f1 value:');
        disp(opti.debug.value(l0_f1));

        disp('l0_s value:');
        disp(opti.debug.value(l0_s));

        disp('l0_f2 value:');
        disp(opti.debug.value(l0_f2));

        disp('T_f1 value:');
        disp(opti.debug.value(T_f1));

        disp('T_s value:');
        disp(opti.debug.value(T_s));

        disp('T_f2 value:');
        disp(opti.debug.value(T_f2));

        disp('Xf1(1,1) value:');
        disp(opti.debug.value(Xf1(1,1)));
    catch err
        disp('Solver failed. Debugging variable values...');
        disp(['T_f1: ', num2str(opti.debug.value(T_f1))]);
        disp(['T_s: ', num2str(opti.debug.value(T_s))]);
        disp(['T_f2: ', num2str(opti.debug.value(T_f2))]);

        disp(['l0_f1: ', num2str(opti.debug.value(l0_f1))]);
        disp(['phi0_f1: ', num2str(opti.debug.value(phi0_f1))]);
        disp(['l0_s: ', num2str(opti.debug.value(l0_s))]);
        disp(['phi0_s: ', num2str(opti.debug.value(phi0_s))]);
        disp(['l0_f2: ', num2str(opti.debug.value(l0_f2))]);
        disp(['phi0_f2: ', num2str(opti.debug.value(phi0_f2))]);

        disp(['Xf1(1,1): ', num2str(opti.debug.value(Xf1(1,1)))]); 
        disp(['Xf1(2,1): ', num2str(opti.debug.value(Xf1(2,1)))]); 
        disp(['Xf1(3,1): ', num2str(opti.debug.value(Xf1(3,1)))]); 
        disp(['Xf1(4,1): ', num2str(opti.debug.value(Xf1(4,1)))]); 
        disp(['Xf1(5,1): ', num2str(opti.debug.value(Xf1(5,1)))]); 
        disp(['Xf1(6,1): ', num2str(opti.debug.value(Xf1(6,1)))]);

        rethrow(err);
    end

    
    % debug minjie
    disp('---- Debugging Variables ----');
    disp(['T_f1: ', num2str(opti.debug.value(T_f1))]);
    disp(['T_s: ', num2str(opti.debug.value(T_s))]);
    disp(['T_f2: ', num2str(opti.debug.value(T_f2))]);

    disp(['l0_f1: ', num2str(opti.debug.value(l0_f1))]);
    disp(['phi0_f1: ', num2str(opti.debug.value(phi0_f1))]);
    disp(['l0_s: ', num2str(opti.debug.value(l0_s))]);
    disp(['phi0_s: ', num2str(opti.debug.value(phi0_s))]);
    disp(['l0_f2: ', num2str(opti.debug.value(l0_f2))]);
    disp(['phi0_f2: ', num2str(opti.debug.value(phi0_f2))]);

    disp(['Xf1(1,1): ', num2str(opti.debug.value(Xf1(1,1)))]);
    disp(['Xf1(2,1): ', num2str(opti.debug.value(Xf1(2,1)))]);
    disp(['Xf1(3,1): ', num2str(opti.debug.value(Xf1(3,1)))]);
    disp(['Xf1(4,1): ', num2str(opti.debug.value(Xf1(4,1)))]);
    disp(['Xf1(5,1): ', num2str(opti.debug.value(Xf1(5,1)))]);
    disp(['Xf1(6,1): ', num2str(opti.debug.value(Xf1(6,1)))]);


    T_f1_sol = sol.value(T_f1);
    T_s_sol  = sol.value(T_s);
    T_f2_sol = sol.value(T_f2);

    disp(['T_f1=', num2str(T_f1_sol), ...
          ', T_s=', num2str(T_s_sol), ...
          ', T_f2=', num2str(T_f2_sol)]);

    Xf1_sol = sol.value(Xf1);
    Xs_sol  = sol.value(Xs);
    Xf2_sol = sol.value(Xf2);


    %Plot
    % time vectors
    tf1_vec = linspace(0, T_f1_sol, N1+1);
    ts_vec  = linspace(T_f1_sol, T_f1_sol+T_s_sol, N2+1);
    tf2_vec = linspace(T_f1_sol+T_s_sol, T_f1_sol+T_s_sol+T_f2_sol, N3+1);

    % plot "foot pos y" vs time
    footY_f1 = zeros(1,N1+1);
    for i=1:(N1+1)
        pp = get_foot_pos(Xf1_sol(:,i), L1, L2);
        footY_f1(i)= pp(2);
    end
    footY_s  = zeros(1,N2+1);
    for i=1:(N2+1)
        pp = get_foot_pos(Xs_sol(:,i), L1, L2);
        footY_s(i)= pp(2);
    end
    footY_f2 = zeros(1,N3+1);
    for i=1:(N3+1)
        pp = get_foot_pos(Xf2_sol(:,i), L1, L2);
        footY_f2(i)= pp(2);
    end

    figure('Name','Foot Y Trajectory'); hold on; grid on;
    plot(tf1_vec, footY_f1, 'r-o','DisplayName','Flight Down');
    plot(ts_vec,  footY_s,  'g-o','DisplayName','Stance');
    plot(tf2_vec, footY_f2, 'b-o','DisplayName','Flight Up');
    legend();
    xlabel('Time (s)'); ylabel('Foot Y (m)');
    title('Foot Y position');

    
    
    % velocity plot
    bodyVelY_f1 = Xf1_sol(4,:);  
    bodyVelY_s  = Xs_sol(4,:);   
    bodyVelY_f2 = Xf2_sol(4,:);  

    figure('Name','Body Y velocity'); 
    hold on; grid on;

    plot(tf1_vec, bodyVelY_f1, 'r-o','DisplayName','Flight Down');
    plot(ts_vec,  bodyVelY_s,  'g-o','DisplayName','Stance');
    plot(tf2_vec, bodyVelY_f2, 'b-o','DisplayName','Flight Up');

    xlabel('Time (s)'); 
    ylabel('Body Velocity in Y (m/s)');
    title('Body Vertical Velocity');
    legend();
    
    % Joint Angles Plot (q1 and q2 over time)
    figure('Name','Joint Angles'); 
    subplot(2,1,1); hold on; grid on;
    plot(tf1_vec, rad2deg(Xf1_sol(2,:)), 'r-o','DisplayName','Flight Down');
    plot(ts_vec,  rad2deg(Xs_sol(2,:)),  'g-o','DisplayName','Stance');
    plot(tf2_vec, rad2deg(Xf2_sol(2,:)), 'b-o','DisplayName','Flight Up');
    xlabel('Time (s)');
    ylabel('q_1 (degrees)');
    title('Joint Angle q_1');
    legend();

    subplot(2,1,2); hold on; grid on;
    plot(tf1_vec, rad2deg(Xf1_sol(3,:)), 'r-o','DisplayName','Flight Down');
    plot(ts_vec,  rad2deg(Xs_sol(3,:)),  'g-o','DisplayName','Stance');
    plot(tf2_vec, rad2deg(Xf2_sol(3,:)), 'b-o','DisplayName','Flight Up');
    xlabel('Time (s)');
    ylabel('q_2 (degrees)');
    title('Joint Angle q_2');
    legend();
    
    %Leg Length l0 over Time 
    figure('Name','Leg Length l0 over Time'); hold on; grid on;
    plot(tf1_vec, l0_f1 * ones(1,N1+1), 'r--','DisplayName','l0 Flight Down');
    plot(ts_vec,  l0_s  * ones(1,N2+1), 'g--','DisplayName','l0 Stance');
    plot(tf2_vec, l0_f2 * ones(1,N3+1), 'b--','DisplayName','l0 Flight Up');
    xlabel('Time (s)'); ylabel('Leg Length l0 (m)');
    title('Leg Length Control');
    legend();

    %Hip Angle phi0 over Time
    figure('Name','Hip Angle phi0 over Time'); hold on; grid on;
    plot(tf1_vec, rad2deg(phi0_f1) * ones(1,N1+1), 'r--','DisplayName','phi0 Flight Down');
    plot(ts_vec,  rad2deg(phi0_s)  * ones(1,N2+1), 'g--','DisplayName','phi0 Stance');
    plot(tf2_vec, rad2deg(phi0_f2) * ones(1,N3+1), 'b--','DisplayName','phi0 Flight Up');
    xlabel('Time (s)'); ylabel('Hip Angle phi0 (degrees)');
    title('Hip Angle Control');
    legend();
    
    %% --- Plot Body Y trajectory ---
    %figure('Name','Body Y Trajectory','Color','white');
    %hold on; grid on;

    % Flight Down:
    %plot(tf1_vec, Xf1_sol(1,:), 'r-o','DisplayName','Flight Down');
    % Stance:
    %plot(ts_vec,  Xs_sol(1,:),  'g-o','DisplayName','Stance');
    % Flight Up:
    %plot(tf2_vec, Xf2_sol(1,:), 'b-o','DisplayName','Flight Up');

    %legend();
    %xlabel('Time (s)');
    %ylabel('Body Y (m)');
    %title('Body Y Position');
end

    