function monopod_1D_three_phase_opt()

    import casadi.*;

    % Duration of the phase
    
    %opti variable T11 ... T3
    T_f1 = 3;  % Flight down
    T_s  = 3;  % Stance
    T_f2 = 3;  % Flight up

    % discrete steps
    N1 = 20; 
    N2 = 20; 
    N3 = 20;

    dt_f1 = T_f1 / N1;
    dt_s  = T_s  / N2;
    dt_f2 = T_f2 / N3;

    g = 9.81; 
    m = 30.0;  % mass
    I1= 0.02; % Joint 1 moment of inertia
    I2= 0.01; % Joint 2 moment of inertia

    opti = Opti();

    % --- Phase1: Flight Down ---
    % state dimension 6: [y, q1, q2, dy, dq1, dq2]
    Xf1 = opti.variable(6, N1+1);
    Uf1 = opti.variable(2, N1);  % torque [u1; u2]

    % --- Phase2: Stance ---
    Xs = opti.variable(6, N2+1);
    Us = opti.variable(2, N2);

    % --- Phase3: Flight Up ---
    Xf2 = opti.variable(6, N3+1);
    Uf2 = opti.variable(2, N3);

   %object function
   % e = Xf2(end) - Xf1(1) soft periodicity condition (as a cost term)
   % 10000 *sumsqr(e)
    obj = sumsqr(Uf1) + sumsqr(Us) + sumsqr(Uf2);

    opti.minimize(obj);

    %Kinetic Constraints: RK4 Multishot 
    
    f_flight = @(x,u) get_dyn_flight_yOnly(x,u,m,g,I1,I2);

    % Phase1: flight down
    for k=1:N1
        xk   = Xf1(:,k);
        uk   = Uf1(:,k);
        xk1  = rk4_step(f_flight, xk, uk, dt_f1);
        opti.subject_to(Xf1(:,k+1) == xk1);
    end

    % Phase3: flight up
    for k=1:N3
        xk   = Xf2(:,k);
        uk   = Uf2(:,k);
        xk1  = rk4_step(f_flight, xk, uk, dt_f2);
        opti.subject_to(Xf2(:,k+1) == xk1);
    end

    % Phase2: stance
    f_stance = @(x,u) get_dyn_stance_yOnly(x,u,m,g,I1,I2);

    for k=1:N2
        xk   = Xs(:,k);
        uk   = Us(:,k);
        xk1  = rk4_step(f_stance, xk, uk, dt_s);
        opti.subject_to(Xs(:,k+1) == xk1);
    end


    % inital conditon: Phase1( flight down )
    opti.subject_to(Xf1(1,1) == 1);   % y=1
    opti.subject_to(Xf1(4,1) == 0);   % dy=0
    opti.subject_to(Xf1(2,1) == deg2rad(-135));   % q1=0
    opti.subject_to(Xf1(3,1) == deg2rad(90));   % q2=0
    opti.subject_to(Xf1(5,1) == 0);   % dq1=0
    opti.subject_to(Xf1(6,1) == 0);   % dq2=0

    % Phase1 end --> Phase2 start
    %Applz impact law     opti.subject_to(Xs(:,1) == get_impact(   Xf1(:,N1+1)) ));

   
    opti.subject_to(Xs(:,1) == Xf1(:,N1+1));
    
    % Touchdown y=0
    %opti.subject_to(Xf1(1,N1+1) == 0);

    % Phase2 end --> Phase3 start
    opti.subject_to(Xf2(:,1) == Xs(:,N2+1));
    %  stance end y=0
    %opti.subject_to(Xs(1,N2+1) == 0);

    % final state: Phase3 end
    opti.subject_to(Xf2(1,N3+1) == 1);  % y=1
    opti.subject_to(Xf2(4,N3+1) == 0);  % dy=0
    opti.subject_to(Xf2(2,N3+1) == deg2rad(-135));  % q1=0
    opti.subject_to(Xf2(3,N3+1) == deg2rad(90));  % q2=0
    opti.subject_to(Xf2(5,N3+1) == 0);  % dq1=0
    opti.subject_to(Xf2(6,N3+1) == 0);  % dq2=0

    % y>=0, Avoid wearing floor
    opti.subject_to(Xf1(1,:) >= 0); 
    opti.subject_to(Xs(1,:)  >= 0);
    opti.subject_to(Xf2(1,:) >= 0);

    % torque limit
    opti.subject_to(-10000 <= Uf1(:) <= 10000); 
    opti.subject_to(-10000 <= Us(:)  <= 10000);
    opti.subject_to(-10000 <= Uf2(:) <= 10000);

    %inital guess
    opti.set_initial(Xf1, 0);
    opti.set_initial(Xs,  0);
    opti.set_initial(Xf2, 0);
    opti.set_initial(Uf1, 0);
    opti.set_initial(Us,  0);
    opti.set_initial(Uf2, 0);

    %solve
    opti.solver('ipopt');
    sol = opti.solve();

    Xf1_sol = sol.value(Xf1);
    Xs_sol  = sol.value(Xs);
    Xf2_sol = sol.value(Xf2);
    Uf1_sol = sol.value(Uf1);
    Us_sol  = sol.value(Us);
    Uf2_sol = sol.value(Uf2);

    %Visualize simple results
    tf1 = linspace(0, T_f1, N1+1);
    ts  = linspace(T_f1, T_f1+T_s, N2+1);
    tf2 = linspace(T_f1+T_s, T_f1+T_s+T_f2, N3+1);

    figure('Name','Trajectory','Color','white');
    subplot(3,1,1);
        hold on; grid on;
        plot(tf1, Xf1_sol(1,:), 'r-o');
        plot(ts,  Xs_sol(1,:),  'g-o');
        plot(tf2, Xf2_sol(1,:), 'b-o');
        xlabel('Time (s)'); ylabel('y');
        title('y Trajectory (Flight-Stance-Flight)');
        legend('FlightDown','Stance','FlightUp');

    subplot(3,1,2);
        hold on; grid on;
        plot(tf1, Xf1_sol(2,:), 'r');
        plot(ts,  Xs_sol(2,:),  'g');
        plot(tf2, Xf2_sol(2,:), 'b');
        xlabel('Time (s)'); ylabel('q1');
        title('q1 angle');

    subplot(3,1,3);
        hold on; grid on;
        plot(tf1, Xf1_sol(3,:), 'r');
        plot(ts,  Xs_sol(3,:),  'g');
        plot(tf2, Xf2_sol(3,:), 'b');
        xlabel('Time (s)'); ylabel('q2');
        title('q2 angle');

    figure('Name','ControlTorques','Color','white');
    subplot(2,1,1);
        hold on; grid on;
        plot(linspace(0,T_f1,N1), Uf1_sol(1,:), 'r');
        plot(linspace(T_f1,T_f1+T_s,N2), Us_sol(1,:), 'g');
        plot(linspace(T_f1+T_s,T_f1+T_s+T_f2,N3), Uf2_sol(1,:), 'b');
        xlabel('Time (s)'); ylabel('u1 (Nm)');
        legend('FlightDown','Stance','FlightUp');

    subplot(2,1,2);
        hold on; grid on;
        plot(linspace(0,T_f1,N1), Uf1_sol(2,:), 'r');
        plot(linspace(T_f1,T_f1+T_s,N2), Us_sol(2,:), 'g');
        plot(linspace(T_f1+T_s,T_f1+T_s+T_f2,N3), Uf2_sol(2,:), 'b');
        xlabel('Time (s)'); ylabel('u2 (Nm)');
end

%RK4 one step 
function x_next = rk4_step(f, x, u, dt)
    k1 = f(x, u);
    k2 = f(x + dt/2*k1, u);
    k3 = f(x + dt/2*k2, u);
    k4 = f(x + dt*k3,   u);
    x_next = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

%dynamic function£ºflight
function dx = get_dyn_flight_yOnly(x,u,m,g,I1,I2)
    % x = [y; q1; q2; dy; dq1; dq2]
    dx = [
        x(4);         % dy
        x(5);         % dq1
        x(6);         % dq2
       -g;            % ddy
        u(1)/I1;      % ddq1
        u(2)/I2       % ddq2
    ];
end

%dynamic function£ºstance 
function dx = get_dyn_stance_yOnly(x,u,m,g,I1,I2)
    % x = [y; q1; q2; dy; dq1; dq2]
    k_leg = 1000;
    y  = x(1);
    dy = x(4);
    F_leg = -k_leg*(y - 0);  
    
    dx = [
        x(4);              % dy
        x(5);              % dq1
        x(6);              % dq2
        F_leg/m - g;       % ddy
        u(1)/I1;           % ddq1
        u(2)/I2            % ddq2
    ];
end
