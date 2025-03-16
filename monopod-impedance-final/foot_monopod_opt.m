function foot_monopod_opt()
    import casadi.*;
    
    opti = casadi.Opti();

    % 3-phase times: T_f1, T_s, T_f2
    T_f1 = opti.variable();  % flight down
    T_s  = opti.variable();  % stance
    T_f2 = opti.variable();  % flight up
    opti.subject_to(T_f1 >= 0);
    opti.subject_to(T_s  >= 0);
    opti.subject_to(T_f2 >= 0);
    opti.subject_to(T_f1 + T_s + T_f2 <= 15);
    

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
    
    % flight down
    K_l = 200;  % example gains
    D_l = 25;
    K_phi = 10;
    D_phi = 10;

    %Decision variables: Xf1, Xs, Xf2; Uf1, Us, Uf2
    % State vector x = [y; q1; q2; dy; dq1; dq2], dimension=6
    Xf1 = opti.variable(6, N1+1);
%     Uf1 = opti.variable(2, N1);   % torque on joint1, joint2

    Xs  = opti.variable(6, N2+1);
%     Us  = opti.variable(2, N2);

    Xf2 = opti.variable(6, N3+1);
%     Uf2 = opti.variable(2, N3);


    % Instead of Uf1, define l0_f1, phi0_f1 minjie
    l0_f1   = opti.variable(1, N1);
    phi0_f1 = opti.variable(1, N1);

    % For stance minjie
    l0_s   = opti.variable(1, N2);
    phi0_s = opti.variable(1, N2);

    % For flight-up minjie
    l0_f2   = opti.variable(1, N3);
    phi0_f2 = opti.variable(1, N3);

    % objective function
    alpha = 1.0;
    e = Xf2(:, N3+1) - Xf1(:,1);   % final - initial
    % For now, only a "periodic" cost:
    obj = alpha * ( sumsqr(e) + 10*sumsqr(Xf2(1, N3+1) - Xf1(1,1)) );
    opti.minimize(obj);
    


    %Build RK4 constraints for each phase
    % flight down
    for k=1:N1
        xk = Xf1(:,k);
%         uk = Uf1(:,k);
%         xk1 = rk4_step(@(xx,uu) get_dyn_flight_foot(xx,uu,m,I1,I2,g), xk, uk, dt_f1);
%         opti.subject_to(Xf1(:,k+1)== xk1);
        l0k   = l0_f1(1,k);
        phi0k = phi0_f1(1,k);
        xk1   = rk4_step(@(xx) get_dyn_flight_imp(xx, l0k, phi0k, ...
                            K_l, D_l, K_phi, D_phi, m, I1, I2, g, L1, L2), ...
                         xk, dt_f1);
        opti.subject_to(Xf1(:,k+1) == xk1);
    end

    % stance
    for k=1:N2
        xk = Xs(:,k);
%         uk = Us(:,k);
%         xk1 = rk4_step(@(xx,uu) get_dyn_stance_foot(xx,uu,m,I1,I2,g), xk, uk, dt_s);
%         opti.subject_to(Xs(:,k+1)== xk1);
        l0k   = l0_s(1, k);         % current leg-length rest
        phi0k = phi0_s(1, k);       % current angle rest

        xk1 = rk4_step(@(xx) get_dyn_stance_imp(xx, l0k, phi0k, ...
                                K_l, D_l, K_phi, D_phi, ...
                                m, I1, I2, g, L1, L2), ...
                       xk, dt_s);

        opti.subject_to(Xs(:, k+1) == xk1);
    end

    % flight up
    for k=1:N3
        xk = Xf2(:,k);
%         uk = Uf2(:,k);
%         xk1 = rk4_step(@(xx,uu) get_dyn_flight_foot(xx,uu,m,I1,I2,g), xk, uk, dt_f2);
%         opti.subject_to(Xf2(:,k+1)== xk1);
        l0k   = l0_f2(1, k);
        phi0k = phi0_f2(1, k);

        xk1 = rk4_step(@(xx) get_dyn_flight_imp(xx, l0k, phi0k, ...
                                K_l, D_l, K_phi, D_phi, ...
                                m, I1, I2, g, L1, L2), ...
                       xk, dt_f2);

        opti.subject_to(Xf2(:, k+1) == xk1);
    end

    % Boundary conditions 
    % initial condition: 
    opti.subject_to(Xf1(1,1) == 1);%y
    %opti.subject_to( 0.9 <= Xf1(1,1) <= 1.1)
    opti.subject_to(Xf1(4,1) == 0);%dy
    opti.subject_to(Xf1(2,1) == deg2rad(-45));%q1
    opti.subject_to(Xf1(3,1) == deg2rad(90));%q2
    opti.subject_to( -0.1 <= Xf1(5,1) <= 0.1)
    opti.subject_to( -0.1 <= Xf1(6,1) <= 0.1)
    %opti.subject_to(Xf1(5,1) >= 0);
    %opti.subject_to(Xf1(6,1) >= 0);
    
    % For example, if your total leg can be about 1.0 m at full extension:
    % Impose a minimum leg length to avoid singularities (e.g. 0.2 m) 
    % and a maximum leg length (e.g. 1.2 m) 
    opti.subject_to( 0.2 <= l0_f1(:) <= 1.0 );
    opti.subject_to( 0.2 <= l0_s(:)  <= 1.0 );
    opti.subject_to( 0.2 <= l0_f2(:) <= 1.0 );

    
    % For angles, let’s assume your leg joint can swing about +/- 90 degrees.
    % That corresponds to +/- pi/2 radians:
    opti.subject_to( -pi/4 <= phi0_f1(:) <= pi/4 );
    opti.subject_to( -pi/4 <= phi0_s(:)  <= pi/4 );
    opti.subject_to( -pi/4 <= phi0_f2(:) <= pi/4 );
    
%     opti.subject_to(Xf2(:, N3+1) == Xf1(:,1));

%     opti.set_initial(l0_f1, 0.5);   
%     opti.set_initial(phi0_f1, 0);
%     opti.set_initial(l0_s,  0.5);
%     opti.set_initial(phi0_s, 0);
%     opti.set_initial(l0_f2,  0.5);
%     opti.set_initial(phi0_f2, 0);??

    %% ==================== Flight Down (l0_f1, phi0_f1) ====================
    % l0_f1, phi0_f1 ????? 1×N1
    % ??? flightDown ??? 0.9 -> 0.7,  phi0 ?? 0

    l0_f1_guess   = zeros(1, N1);
    phi0_f1_guess = zeros(1, N1);

    for k=1:N1
        alpha = (k-1)/(N1-1);  % alpha in [0..1]
        % ? l0_f1 ? 0.9 ????? 0.7
        l0_f1_guess(k)   = 0.9 - 0.2*alpha;

        % phi0_f1_guess ?? 0, ??????
        phi0_f1_guess(k) = 0.0;  
    end

    opti.set_initial(l0_f1,   l0_f1_guess);
    opti.set_initial(phi0_f1, phi0_f1_guess);

    %% ==================== Stance (l0_s, phi0_s) ============================
    % l0_s, phi0_s ????? 1×N2
    % ??? stance ?? l0_s ????? 0.7, phi0_s ???? 0

    l0_s_guess   = 0.7 * ones(1, N2);
    phi0_s_guess = 0.0 * ones(1, N2);

    opti.set_initial(l0_s,   l0_s_guess);
    opti.set_initial(phi0_s, phi0_s_guess);

    %% ==================== Flight Up (l0_f2, phi0_f2) ======================
    % l0_f2, phi0_f2 ????? 1×N3
    % ??? flightUp ??? 0.7 -> 0.9, phi0_f2?? 0

    l0_f2_guess   = zeros(1, N3);
    phi0_f2_guess = zeros(1, N3);

    for k=1:N3
        alpha = (k-1)/(N3-1);
        l0_f2_guess(k)   = 0.65 + 0.3*alpha;  
        phi0_f2_guess(k) = 0.0;
    end

    opti.set_initial(l0_f2,   l0_f2_guess);
    opti.set_initial(phi0_f2, phi0_f2_guess);



    % Phase1 -> Phase2 continuity at foot collision (foot y=0)
    pFoot_f1_end = get_foot_pos(Xf1(:,N1+1), L1, L2);
    opti.subject_to( pFoot_f1_end(2) >= -0.05 );  % touch ground
    opti.subject_to( pFoot_f1_end(2) <= 0.05 );

    % continuity in states
    opti.subject_to( Xs(:,1) == Xf1(:,N1+1));

    % Phase2 -> Phase3 continuity => foot y>0 => leaving ground
    opti.subject_to( Xf2(:,1) == Xs(:,N2+1));
    
    %path constraints: y >= 0 => body doesn't go below ground?
    for k=1:(N1+1)
       footPos_k = get_foot_pos(Xf1(:,k), L1, L2);
       opti.subject_to( footPos_k(2) >= 0 );  
       %opti.subject_to( Xf1(1,k) >= 0 );
    end
    for k=1:(N3+1)
       footPos_k = get_foot_pos(Xf2(:,k), L1, L2);
       opti.subject_to( footPos_k(2) >= -0.05 );
       %opti.subject_to( Xf2(1,k) >= 0 );
    end
    % In stance, footPos=0 => keep foot locked? 
    for k=1:(N2+1)
       footPos_k = get_foot_pos(Xs(:,k), L1, L2);
       opti.subject_to( -0.05 <= footPos_k(2) <= 0.05 );
       %opti.subject_to( Xs(1,k) >= 0 );
    end

    % torque limit
    maxTrq = 5000;
%     opti.subject_to( -maxTrq <= Uf1(:) <= maxTrq );
%     opti.subject_to( -maxTrq <= Us(:)  <= maxTrq );
%     opti.subject_to( -maxTrq <= Uf2(:) <= maxTrq );

    %initial guess
%     opti.set_initial(Xf1, 0);
%     opti.set_initial(Xs,  0);
%     opti.set_initial(Xf2, 0);??
%     opti.set_initial(Uf1, 0);
%     opti.set_initial(Us,  0);
%     opti.set_initial(Uf2, 0);
    % ========= 1) Flight Down (Xf1) from y=1 => y=0.3 ======================
    Xf1_guess = zeros(6, N1+1);  % 6? × (N1+1)?

    for k = 1:(N1+1)
        alpha = (k-1)/N1;              % alpha?0?1
        % --- y?1????0.3 ---
        y_val   = 1.0 - 0.7 * alpha;   % 1 -> 0.3
        % --- q1,q2?????(-45°, 90°) ---
        q1_val  = deg2rad(-45);
        q2_val  = deg2rad(90);
        % --- dy??????, ?????? ---
        dy_val  = -0.5;  % ????????0?-0.5
        % --- dq1,dq2??0 ---
        dq1_val = 0.0;
        dq2_val = 0.0;

        Xf1_guess(:, k) = [y_val; q1_val; q2_val; dy_val; dq1_val; dq2_val];
    end

    % ?Xf1_guess??Xf1??
    opti.set_initial(Xf1, Xf1_guess);


    % ========= 2) Stance (Xs) from y=0.3 ???? =========================
    Xs_guess = zeros(6, N2+1);

    for k = 1:(N2+1)
        alpha = (k-1)/N2;              % 0..1
        % --- y???0.3 ---
        y_val   = 0.3;
        % --- ??q1?-45°???-30°, q2?90°?60°?? ---
        q1_val  = deg2rad(-45) + alpha*(deg2rad(-30)-deg2rad(-45));
        q2_val  = deg2rad(90) + alpha*(deg2rad(60)-deg2rad(90));
        % --- dy???0, ??????? ---
        dy_val  = 0.0;
        % --- dq1,dq2???0 ---
        dq1_val = 0.0;
        dq2_val = 0.0;

        Xs_guess(:, k) = [y_val; q1_val; q2_val; dy_val; dq1_val; dq2_val];
    end

    opti.set_initial(Xs, Xs_guess);


    % ========= 3) Flight Up (Xf2) from y=0.3 => y=1.0 ======================
    Xf2_guess = zeros(6, N3+1);

    for k = 1:(N3+1)
        alpha = (k-1)/N3;               % 0..1
        % --- y?0.3?????1.0 ---
        y_val   = 0.3 + 0.7 * alpha;    % 0.3 -> 1.0
        % --- q1,q2???? stance???????,?????? ---
        q1_stanceEnd = deg2rad(-30);    % ??stanc???
        q2_stanceEnd = deg2rad(60);
        q1_val  = q1_stanceEnd; 
        q2_val  = q2_stanceEnd; 
        % --- dy???????, ???? ---
        dy_val  = 1.0;   
        % --- dq1,dq2??0 ---
        dq1_val = 0.0;
        dq2_val = 0.0;

        Xf2_guess(:, k) = [y_val; q1_val; q2_val; dy_val; dq1_val; dq2_val];
    end

    opti.set_initial(Xf2, Xf2_guess);

    opti.set_initial(T_f1, 1);
    opti.set_initial(T_s,  0.2);
    opti.set_initial(T_f2, 1);
    
    %final state equals initial state 
    %opti.subject_to(Xf2(1, end) == 1);
    %opti.subject_to(Xf2(4, end) == 0);
    %opti.subject_to(Xf2(2, end) == deg2rad(-135));
    %opti.subject_to(Xf2(3, end) == deg2rad(90));
    %opti.subject_to(Xf2(5, end) == 0);
    %opti.subject_to(Xf2(6, end) == 0);

    %Solve
    opts = struct;
    opts.ipopt.max_iter = 3000;
    opts.ipopt.acceptable_tol = 1e-1;  % ?????1e-6
    opts.ipopt.acceptable_constr_viol_tol = 1e-1;
    opts.ipopt.acceptable_dual_inf_tol = 1e+2; 
    opts.ipopt.dual_inf_tol            = 1e+2;
%     opts.ipopt.constr_viol_tol = 1e-1; 
%     opts.ipopt.dual_inf_tol = 1e-1; 
%     opts.ipopt.acceptable_tol = 1e-1; 
%     opts.ipopt.acceptable_constr_viol_tol = 2e-1;
    opti.solver('ipopt', opts);
    sol = opti.solve();

    T_f1_sol = sol.value(T_f1);
    T_s_sol  = sol.value(T_s);
    T_f2_sol = sol.value(T_f2);

    disp(['T_f1=', num2str(T_f1_sol), ...
          ', T_s=', num2str(T_s_sol), ...
          ', T_f2=', num2str(T_f2_sol)]);

    Xf1_sol = sol.value(Xf1);
    Xs_sol  = sol.value(Xs);
    Xf2_sol = sol.value(Xf2);
%     Uf1_sol = sol.value(Uf1);
%     Us_sol  = sol.value(Us);
%     Uf2_sol = sol.value(Uf2);

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

    % torque plot
    t_f1 = linspace(0, T_f1_sol, N1);
    t_s  = linspace(T_f1_sol, T_f1_sol+T_s_sol, N2);
    t_f2 = linspace(T_f1_sol+T_s_sol, T_f1_sol+T_s_sol+T_f2_sol, N3);
% 
%     figure('Name','ControlTorques');
%     subplot(2,1,1); hold on; grid on;
%     plot(t_f1, Uf1_sol(1,:), 'r');
%     plot(t_s,  Us_sol(1,:),  'g');
%     plot(t_f2, Uf2_sol(1,:), 'b');
%     ylabel('u_1 (Nm)');
%     legend('FlightDown','Stance','FlightUp');
% 
%     subplot(2,1,2); hold on; grid on;
%     plot(t_f1, Uf1_sol(2,:), 'r');
%     plot(t_s,  Us_sol(2,:),  'g');
%     plot(t_f2, Uf2_sol(2,:), 'b');
%     xlabel('Time (s)'); ylabel('u_2 (Nm)');
%     legend('FlightDown','Stance','FlightUp');
    
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

    