%% Simulate monopod dynamics
% Init
clear; clc; close all;

% Toggle simulation behavior
hopping_control_mode = 2;  % 0 = No Hopping Control, 1 = Raibert, 2 = Feedback, 3 = Impulse
impedance_control_mode = 3;  % 1 = Joint Impedance, 2 = Cartesian Impedance, 3 = Polar Impedance 
export_graphs = false;

% Set desired hopping heights (foot) and corresponding times
H_d = [0.3 0.6 0.1];
t_d = [0 5 10];

d_counter = 1;

% State vector: x = [y; q1; q2; dy; dq1; dq2]
% Control vector: u = [u1; u2]

x_0 = [0.23, deg2rad(-135), deg2rad(90), 0, 0, 0];

x_full = x_0;
t_full = 0;

t_max = 15;
t_step = 0.025;
t_span = 0:t_step:t_max;

l_rest = get_l_rest();

% Event functions
opts_flight1 = odeset('Events', @guard_high);
opts_flight2 = odeset('Events', @guard_TD);
opts_stance1 = odeset('Events', @guard_low);

switch hopping_control_mode
    case 0  % No Hopping Control
        l_d = l_rest;
    case 2  % Feedback
        r_F_pol = get_r_F_pol([1, deg2rad(-135), deg2rad(90), 0, 0, 0]');
        l_d = r_F_pol(1);
        x_high_prev = x_0;
    case 3  % Impulse
        l_d = l_rest;
        x_high_prev = x_0;
end

while true    
    % Flight phase 2: From high-point to touchdown
    [t_out_f2, x_out_f2, te_f2, xe_f2, ie_f2] = ode45(@(t, x) get_dyn_flight(t, x, impedance_controller(x, [l_rest; -pi/2], impedance_control_mode)), t_span, x_full(end, :), opts_flight2);
    x_full = [x_full; x_out_f2(2:end, :)];  % remove first entry of x_out since it is already in last entry of x_full
    t_full = [t_full; t_out_f2(2:end)];

    t_span = t_full(end):t_step:t_max;

    if ie_f2
        fprintf("t = %.2fs: Touchdown", t_full(end))
    end
    
    if t_full(end) >= t_max - t_step
        break
    end

    % Impact
    x_full(end, 4:6) = get_impact(xe_f2')';
    
    if hopping_control_mode == 1
        delta_l_loss = get_delta_l_loss(xe_f2');
        l_d = l_rest + delta_l_loss;
    end
    
    fprintf(" & Impact\n")
    
    if t_full(end) >= t_max - t_step
        break
    end

    % Stance phase 1: From touchdown to low-point
    [t_out_s1, x_out_s1, te_s1, xe_s1, ie_s1] = ode45(@(t, x) get_dyn_stance(t, x, impedance_controller(x, [l_d; -pi/2], impedance_control_mode)), t_span, x_full(end, :), opts_stance1);
    x_full = [x_full; x_out_s1(2:end, :)];
    t_full = [t_full; t_out_s1(2:end)];

    t_span = t_full(end):t_step:t_max;

    if ie_s1
        fprintf("t = %.2fs: Low-point reached\n", t_full(end))
    end

    if t_full(end) >= t_max - t_step
        break
    end
    
    % Set desired hopping height corresponding to current time
    if t_full(end) >= t_d(d_counter)
        H = H_d(d_counter);
        d_counter = min(d_counter + 1, length(t_d)); 
    end
    
    switch hopping_control_mode           
        case 0  % No Hopping Control
            l_d = l_rest;

        case 1  % Raibert
            delta_l_height = get_delta_l_height(xe_s1', H);
            l_d = l_rest + delta_l_loss + delta_l_height;

        case 2  % Feedback
            delta_l_height = get_delta_l(x_high_prev', H);
            l_d = l_d + delta_l_height;

        case 3  % Impulse
            delta_l_height = impulse_hopping_height_controller(H, x_high_prev(1));
            l_d = l_d + delta_l_height;
    end

    % Stance phase 2: From low-point to lift-off
    opts_stance2 = odeset('Events', @(t,x) guard_LO(t, x, [l_d; -pi/2], impedance_control_mode));
    [t_out_s2, x_out_s2, te_s2, xe_s2, ie_s2] = ode45(@(t, x) get_dyn_stance(t, x, impedance_controller(x, [l_d; -pi/2], impedance_control_mode)), t_span, x_full(end, :), opts_stance2);
    x_full = [x_full; x_out_s2(2:end, :)];
    t_full = [t_full; t_out_s2(2:end)];

    t_span = t_full(end):t_step:t_max;

    if ie_s2
        fprintf("t = %.2fs: Lift-Off\n", t_full(end))
    end

    if t_full(end) >= t_max - t_step
        break
    end
    
    % Flight phase 1: From lift-off to high-point
    [t_out_f1, x_out_f1, te_f1, xe_f1, ie_f1] = ode45(@(t, x) get_dyn_flight(t, x, impedance_controller(x, [l_rest; -pi/2], impedance_control_mode)), t_span, x_full(end, :), opts_flight1);
    x_high_prev = xe_f1;
    x_full = [x_full; x_out_f1(2:end, :)];
    t_full = [t_full; t_out_f1(2:end)];

    t_span = t_full(end):t_step:t_max;

    if ie_f1
        fprintf("t = %.2fs: High-point reached\n", t_full(end))
    end
    
    if t_full(end) >= t_max - t_step
        break
    end
end

% Plot trajectory
monoplot(t_full, x_full, t_d, H_d, impedance_control_mode, hopping_control_mode, export_graphs)
% plot_foot(t_full, x_full, t_d, H_d, hopping_control_mode, export_graphs)
% plot_body(t_full, x_full, t_d, H_d, hopping_control_mode, export_graphs)
