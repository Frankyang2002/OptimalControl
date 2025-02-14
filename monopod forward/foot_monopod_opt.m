function foot_monopod_opt()
    import casadi.*;
    opti = casadi.Opti();
    % f1 = falling, s = stance, f2 = jump up fly

    % 3-phase times: T_f1, T_s, T_f2
    T_f1 = opti.variable();  % flight down
    T_s  = opti.variable();  % stance
    T_f2 = opti.variable();  % flight up

    % Decrease time space to search smaller space
    opti.subject_to(0.1 <= T_f1 <= 2.5); 
    opti.subject_to(0.1 <= T_s <= 2.0);  
    opti.subject_to(0.1 <= T_f2 <= 2.5); 
    opti.subject_to(T_f1 + T_s + T_f2 <= 5.0); 

    % number of intervals each phase
    N1 = 20; N2 = 20; N3 = 20;
    dt_f1 = T_f1 / N1;
    dt_s  = T_s  / N2;
    dt_f2 = T_f2 / N3;

    %Model parameters 
    m     = 80;   % body mass
    I1    = 0.01; % inertia of joint1
    I2    = 0.01; % inertia of joint2
    g     = 9.81; 
    L1    = 0.5;  % link1 length
    L2    = 0.5;  % link2 length
    LRest = 0.70;   % around 1/sqrt(2) from base arrangement of 45 90 degrees
    k_leg = 5500;

    %Decision variables: Xf1, Xs, Xf2; Uf1, Us, Uf2
    % State vector x = [x; y; q1; q2; dx; dy; dq1; dq2], dimension=8
    Xf1 = opti.variable(8, N1+1);
    Uf1 = opti.variable(2, N1);   % torque on joint1, joint2

    Xs  = opti.variable(8, N2+1);
    Us  = opti.variable(2, N2);

    Xf2 = opti.variable(8, N3+1);
    Uf2 = opti.variable(2, N3);

    %Objective function: sum of torque^2 + soft periodic
    alpha = 1e4;
    e = Xf2(:, N3+1) - Xf1(:,1);  % final state - initial state, Punish state errors
    obj = sumsqr(Uf1) + sumsqr(Us) + sumsqr(Uf2) + alpha*sumsqr(e);
    opti.minimize(obj);

    %Build RK4 constraints for each phase
    % flight down
    for k=1:N1
        xk = Xf1(:,k);
        uk = Uf1(:,k);
        xk1 = rk4_step(@(xx,uu) get_dyn_flight_foot(xx,uu,m,I1,I2,g), xk, uk, dt_f1);
        opti.subject_to(Xf1(:,k+1)== xk1);
    end
 
    % stance
    for k=1:N2
        xk = Xs(:,k);
        uk = Us(:,k);
        xk1 = rk4_step(@(xx,uu) get_dyn_stance_foot(xx,uu,m,I1,I2,g,L1,L2,LRest,k_leg), xk, uk, dt_s);
        opti.subject_to(Xs(:,k+1)== xk1);
    end

    % flight up
    for k=1:N3
        xk = Xf2(:,k);
        uk = Uf2(:,k);
        xk1 = rk4_step(@(xx,uu) get_dyn_flight_foot(xx,uu,m,I1,I2,g), xk, uk, dt_f2);
        opti.subject_to(Xf2(:,k+1)== xk1);
    end

    %% --- Initial Conditions ---
    % x = [x; y; q1; q2; dx; dy; dq1; dq2]
    opti.subject_to(Xf1(1,1) == 0);%x starts at 0
    opti.subject_to(Xf1(2,1) == 1);%y, starts at 1m
    
    opti.subject_to(Xf1(3,1) == deg2rad(90));%q1 
    opti.subject_to(Xf1(4,1) == deg2rad(270));%q2 

    opti.subject_to(Xf1(5,1) >= 0.5);%dx, Initial x vel
    opti.subject_to(Xf1(6,1) == 0);%dy, no initial vertical vel

    % limit angle change at start
    opti.subject_to(-0.1 <= Xf1(7:8,1) <= 0.1)




    %% --- Guard Function ---
    % Phase1 -> Phase2 continuity at foot collision (foot y=0)
    pFoot_f1_end = get_foot_pos(Xf1(:,N1+1), L1, L2);
    opti.subject_to( pFoot_f1_end(2) == 0 );  % touch ground

    % continuity in states
    opti.subject_to( Xs(:,1) == Xf1(:,N1+1));

    % Phase2 -> Phase3 continuity => foot y>0 => leaving ground 
    opti.subject_to( Xf2(:,1) == Xs(:,N2+1));



    %% --- Path Constraints ---
    %path constraints: y >= ymin => body doesn't rouch too far down
    ymin = 0.1;
    opti.subject_to(Xf1(2, :) >= ymin);
    opti.subject_to(Xs(2, :) >= ymin);
    opti.subject_to(Xf2(2, :) >= ymin);

    % max distance
    opti.subject_to(Xf2(1, N3+1) <= 20);

    % Foot does not go underground
    for k=1:(N1+1)
       footPos_k = get_foot_pos(Xf1(:,k), L1, L2);
       opti.subject_to( footPos_k(2) >= 0 );  
    end

    % Stance foot already checked with the dynamics file 
    for k=1:(N1+1)
       footPos_k = get_foot_pos(Xs(:,k), L1, L2);
       opti.subject_to( footPos_k(2) == 0 );  
    end

    % Flight 2, make sure its above ground
    for k=1:(N3+1)
       footPos_k = get_foot_pos(Xf2(:,k), L1, L2);
       opti.subject_to( footPos_k(2) >= 0 );
    end
    



    % Make sure it is moving
    % x = [x; y; q1; q2; dx; dy; dq1; dq2]
    opti.subject_to(0.5 <= Xf1(5,:) <= 5);
    opti.subject_to(0.5 <= Xs(5,:) <= 5);
    opti.subject_to(0.5 <= Xf2(5,:) <= 5);

    % Angle limit
    % leg doesnt stick up into the body
    opti.subject_to(0 <= Xf1(3,:) <= pi)
    opti.subject_to(0 <= Xs(3,:) <= pi)
    opti.subject_to(0 <= Xf2(3,:) <= pi)

    % 2nd leg can fold in or be straight like a real leg
    opti.subject_to(pi <= Xf1(4,:) <= 2*pi)
    opti.subject_to(pi <= Xs(4,:) <= 2*pi)
    opti.subject_to(pi <= Xf2(4,:) <= 2*pi)


    % torque limit
    maxTrq = 10000;
    opti.subject_to( -maxTrq <= Uf1(:) <= maxTrq );
    opti.subject_to( -maxTrq <= Us(:)  <= maxTrq );
    opti.subject_to( -maxTrq <= Uf2(:) <= maxTrq );

    % Positions
    opti.set_initial(Xf1(1,:), 0);  % y
    opti.set_initial(Xf1(2,:), 1.5);  % y
    
    % Angles
    opti.set_initial(Xf1(3,:), deg2rad(45));             % q1
    opti.set_initial(Xf1(4,:), deg2rad(240));            % q2 straight leg
    
    % Velocities
    opti.set_initial(Xf1(5,:), 1.0);                     % dx
    opti.set_initial(Xf1(6,:), 0);                       % dy
    opti.set_initial(Xf1(7,:), 0);                       % dq1 ≈ 0 rad/s
    opti.set_initial(Xf1(8,:), 0);                       % dq2 ≈ 0 rad/s



    % time
    opti.set_initial(T_f1, 0.5);
    opti.set_initial(T_s,  0.2);
    opti.set_initial(T_f2, 0.5);
    



    %% ---  Solve ---
    opti.solver('ipopt');
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
    Uf1_sol = sol.value(Uf1);
    Us_sol  = sol.value(Us);
    Uf2_sol = sol.value(Uf2);







    %% --- Plots ---
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
    
    %{
    % torque plot
    t_f1 = linspace(0, T_f1_sol, N1);
    t_s  = linspace(T_f1_sol, T_f1_sol+T_s_sol, N2);
    t_f2 = linspace(T_f1_sol+T_s_sol, T_f1_sol+T_s_sol+T_f2_sol, N3);

    
    
    figure('Name','ControlTorques');
    subplot(2,1,1); hold on; grid on;
    plot(t_f1, Uf1_sol(1,:), 'r');
    plot(t_s,  Us_sol(1,:),  'g');
    plot(t_f2, Uf2_sol(1,:), 'b');
    ylabel('u_1 (Nm)');
    legend('FlightDown','Stance','FlightUp');

    subplot(2,1,2); hold on; grid on;
    plot(t_f1, Uf1_sol(2,:), 'r');
    plot(t_s,  Us_sol(2,:),  'g');
    plot(t_f2, Uf2_sol(2,:), 'b');
    xlabel('Time (s)'); ylabel('u_2 (Nm)');
    legend('FlightDown','Stance','FlightUp');


    % x = [x; y; q1; q2; dx; dy; dq1; dq2]
    % velocity plot
    bodyVelY_f1 = Xf1_sol(6,:);  
    bodyVelY_s  = Xs_sol(6,:);   
    bodyVelY_f2 = Xf2_sol(6,:);  

    figure('Name','Body Y velocity'); 
    hold on; grid on;

    plot(tf1_vec, bodyVelY_f1, 'r-o','DisplayName','Flight Down');
    plot(ts_vec,  bodyVelY_s,  'g-o','DisplayName','Stance');
    plot(tf2_vec, bodyVelY_f2, 'b-o','DisplayName','Flight Up');

    % X
    xlabel('Time (s)'); 
    ylabel('Body Velocity in X (m/s)');
    title('Body Vertical Velocity');
    legend();

    bodyVelX_f1 = Xf1_sol(5,:);  
    bodyVelX_s  = Xs_sol(5,:);   
    bodyVelX_f2 = Xf2_sol(5,:);  

    figure('Name','Body X velocity'); 
    hold on; grid on;

    plot(tf1_vec, bodyVelX_f1, 'r-o','DisplayName','Flight Down');
    plot(ts_vec,  bodyVelX_s,  'g-o','DisplayName','Stance');
    plot(tf2_vec, bodyVelX_f2, 'b-o','DisplayName','Flight Up');

    xlabel('Time (s)'); 
    ylabel('Body Velocity in X (m/s)');
    title('Body Horizontal Velocity');
    legend();
    
    
    % Joint Angles Plot (q1 and q2 over time)
    figure('Name','Joint Angles'); 
    subplot(2,1,1); hold on; grid on;
    plot(tf1_vec, rad2deg(Xf1_sol(3,:)), 'r-o','DisplayName','Flight Down');
    plot(ts_vec,  rad2deg(Xs_sol(3,:)),  'g-o','DisplayName','Stance');
    plot(tf2_vec, rad2deg(Xf2_sol(3,:)), 'b-o','DisplayName','Flight Up');
    xlabel('Time (s)');
    ylabel('q_1 (degrees)');
    title('Joint Angle q_1');
    legend();

    subplot(2,1,2); hold on; grid on;
    plot(tf1_vec, rad2deg(Xf1_sol(4,:)), 'r-o','DisplayName','Flight Down');
    plot(ts_vec,  rad2deg(Xs_sol(4,:)),  'g-o','DisplayName','Stance');
    plot(tf2_vec, rad2deg(Xf2_sol(4,:)), 'b-o','DisplayName','Flight Up');
    xlabel('Time (s)');
    ylabel('q_2 (degrees)');
    title('Joint Angle q_2');
    legend();
    %}
   
    %% --- Plot Body xy trajectory ---
    figure('Name','Body Y Trajectory','Color','white');
    hold on; grid on;

    % Flight Down:
    plot(Xf1_sol(1,:), Xf1_sol(2,:), 'r-o','DisplayName','Flight Down');
    % Stance:
    plot(Xs_sol(1,:),  Xs_sol(2,:),  'g-o','DisplayName','Stance');
    % Flight Up:
    plot(Xf2_sol(1,:), Xf2_sol(2,:), 'b-o','DisplayName','Flight Up');

    legend();
    xlabel('Body X (m)');
    ylabel('Body Y (m)');
    title('Body XY Position');
end

    