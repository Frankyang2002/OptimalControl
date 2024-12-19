function [value, isterminal, direction] = guard_TD(t, x)
    % Touchdown occurs when the foot y-coordinate equals zero
    value      = get_c(x);
    isterminal = 1;
    direction  = -1;
end
