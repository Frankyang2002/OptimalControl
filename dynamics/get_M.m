function M = get_M(in1)
%GET_M
%    M = GET_M(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    23-Nov-2022 15:51:41

q1 = in1(2,:);
q2 = in1(3,:);
t2 = cos(q1);
t3 = cos(q2);
t4 = q1+q2;
t5 = cos(t4);
t6 = t2.*1.688e-2;
t7 = t3.*3.968e-4;
t8 = t5.*2.48e-3;
t9 = t7+1.984e-4;
t10 = t6+t8;
M = reshape([3.29e+2./1.0e+3,t10,t8,t10,t3.*7.936e-4+1.9456e-3,t9,t8,t9,6.045916666666667e-4],[3,3]);
