function J_pq = get_J_pq(x)
    q1 = x(2);
    q2 = x(3);
    
    % Link lengths
    L1 = 0.5;
    L2 = 0.5;
    
    % Compute partial derivatives
    dl_dq1 = L1 * cos(q1) + L2 * cos(q1 + q2);
    dl_dq2 = L2 * cos(q1 + q2);
    
    epsilon = 1e-6;
    dphi_dq1 = -sin(q1) / (cos(q1) + cos(q1 + q2) + epsilon);
    dphi_dq2 = -sin(q1 + q2) / (cos(q1) + cos(q1 + q2) + epsilon);

    
    % Construct Jacobian matrix
    J_pq = [dl_dq1, dl_dq2;
            dphi_dq1, dphi_dq2];

    % Debugging: Display size
    disp('Checking J_pq size:');
    disp(size(J_pq));
end
