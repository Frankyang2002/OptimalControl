function out1 = get_J_v_LL(in1)
%GET_J_V_LL
%    OUT1 = GET_J_V_LL(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    23-Nov-2022 15:51:41

q1 = in1(2,:);
q2 = in1(3,:);
t2 = q1+q2;
t3 = cos(t2);
t4 = sin(t2);
t5 = t3.*(2.0./2.5e+1);
t6 = t4.*(2.0./2.5e+1);
t7 = -t6;
out1 = reshape([0.0,1.0,0.0,t7-sin(q1).*(4.0./2.5e+1),t5+cos(q1).*(4.0./2.5e+1),0.0,t7,t5,0.0],[3,3]);