function [value, isterminal, direction] = guard_high(t, x)    
    % High point is reached when y velocity changes from positiv to negative
    dy = x(4);
    
    value      = dy;
    isterminal = 1;
    direction  = -1;
end
