function [value, isterminal, direction] = guard_LO(t, x, x_set, control_mode)
    % Lift-off occurs when constraint forces in y-direction equal zero
    
    u = impedance_controller(x, x_set, control_mode);
    lambda = get_lambda(x, u);
    F_y = lambda(2);
    
    value      = F_y;
    isterminal = 1;
    direction  = -1;
end
