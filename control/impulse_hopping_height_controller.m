%% Impulse Hopping Height Controller
function delta_l = impulse_hopping_height_controller(y_set, y_apex)
    
    % Controller Parameters
    %t_max = 20 * t_step;
    k = 10;

    %t = 0:t_step:t_max;

    % Calculate the height error
    e = y_set - y_apex;
    
    % Set parameters for impulse function
    %mu = t_max/2;
    %sig = 5 * t_step;
    %amplitude = k * e;
    %parameters = [mu, sig, amplitude];

    delta_l = k * e;

end
