function J_pq = get_J_pq(in1)
%GET_J_PQ
%    J_PQ = GET_J_PQ(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    23-Nov-2022 15:51:40

q2 = in1(3,:);
J_pq = reshape([0.0,1.0,sqrt(2.0).*sin(q2).*1.0./sqrt(cos(q2)+1.0).*(-2.0./2.5e+1),1.0./2.0],[2,2]);
