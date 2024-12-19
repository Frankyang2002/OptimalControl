function monoplot(t, x, t_d, H_d, impedance_control_mode, hopping_control_mode, export)
    load('monopod_parameters', 'l_B', 'l_UL', 'l_LL', 'w_B', 'w_UL', 'w_LL')
    
    t_step = t(end)/length(t);
    
    figure
    xlim([-0.7 0.7])
    ylim([-0.1 1.3])
    line([-0.7 0.7], [0 0], 'color', 'b')
    l = line(0, 0);
    
    color = [0.8 0.8 0.8];
    p_B = patch(zeros(4, 1), zeros(4, 1), color);
    p_UL = patch(zeros(4, 1), zeros(4, 1), color);
    p_LL = patch(zeros(4, 1), zeros(4, 1), color);
    
    counter = 1;
    
    if export
        switch hopping_control_mode
            case 0  % No Hopping Control
                switch impedance_control_mode
                    case 1
                        file_name = 'monoplot_impact_Joint_Impedance.gif';
                    case 2
                        file_name = 'monoplot_impact_Cartesian_Impedance.gif';
                    case 3
                        file_name = 'monoplot_impact_Polar_Impedance.gif';
                end
            case 1
                file_name = 'monoplot_Raibert.gif';
            case 2
                file_name = 'monoplot_Feedback.gif';
            case 3
                file_name = 'monoplot_Impulse.gif';
        end
    
    
        if exist('img', 'dir')
            path_dest = append(pwd, filesep, 'img', filesep, file_name);
        else
            path_dest = append(pwd, filesep);
        end
    end
    
    for i = 1:length(t)
        tic
        
        q = x(i, 1:3)';
        
        p_B = transform_patch(p_B, l_B, w_B, get_W_X_B(q));
        p_UL = transform_patch(p_UL, l_UL, w_UL, get_W_X_UL(q));
        p_LL = transform_patch(p_LL, l_LL, w_LL, get_W_X_LL(q));
        
        if t(i) >= t_d(counter)
            H = H_d(counter);
            counter = min(counter + 1, length(t_d)); 
        end
        
        if hopping_control_mode
            delete(l)
            l = line([-0.7 0.7], [H H], 'color', 'r');
        end
        
        title(sprintf("t = %.2fs", t(i)))
        
        if export && i == 1
            gif(path_dest, 'DelayTime', t_step, 'overwrite', true)
        elseif export
            gif
        end
        
        t_elapsed = toc;
        pause(max(t_step - t_elapsed, 0))  % adjust pause such that plot updates in (near) real-time
    end
end


function p_trans = transform_patch(p, length, width, transform)
    % initialize vertices around world origin and transform to specified 
    % pose using W_X_*, update existing patch with transformed vertices
    vertices_homo = [-length/2, -length/2, length/2, length/2;
                      -width/2,   width/2,  width/2, -width/2;
                             0,         0,        0,        0;
                             1,         1,        1,        1];
    
    vertices_trans = transform * vertices_homo;
    p_trans = p;
    p_trans.Vertices = vertices_trans(1:2, :)';
end
