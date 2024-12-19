function n = get_n(in1,in2)
%GET_N
%    N = GET_N(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    23-Nov-2022 15:51:41

dq1 = in2(2,:);
dq2 = in2(3,:);
dy = in2(1,:);
q1 = in1(2,:);
q2 = in1(3,:);
t2 = sin(q2);
t3 = q1+q2;
t4 = cos(t3);
t5 = sin(t3);
t6 = t4.*2.43288e-2;
n = [dq2.*dy.*t5.*(-4.96e-3)-(dq1.*dy.*(t5.*6.2e+1+sin(q1).*2.11e+2))./1.25e+4+3.22749;t6+cos(q1).*1.655928e-1-dq1.^2.*t2.*3.968e-4-dq1.*dq2.*t2.*7.936e-4;t6+dq1.*dq2.*t2.*3.968e-4];