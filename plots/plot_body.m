function plot_body(t, x, t_d, H_d, control_mode, export)
    y = [];
    
    for i = 1:length(t)
        q = x(i, 1:3)';
        W_X_B = get_W_X_B(q);
        
        y = [y W_X_B(2, 4)];
    end
    
    f = figure;
    hold on
    xlim([0 t(end)])
    ylim([-0.1 1.3])
    plot(t, y)
    
    % Plot desired values
    t_plot = [t_d(1)];
    t_plot = [t_plot repelem(t_d(2:end), 2)];
    t_plot = [t_plot t(end)];
    H_plot = repelem(H_d, 2);

    plot(t_plot, H_plot, 'color', 'r')
    
    if export
        switch control_mode
            case 1
                file_name = 'y_body_Raibert.png';
            case 2
                file_name = 'y_body_Feedback.png';
            case 3
                file_name = 'y_body_Impulse.png';
        end

        if exist('img', 'dir')
            path_dest = append(pwd, filesep, 'img', filesep, file_name);
        else
            path_dest = append(pwd, filesep);
        end

        exportgraphics(f, path_dest)
    end
end
