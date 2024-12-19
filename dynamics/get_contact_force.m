function lambda = get_contact_force(q, dq)


k = 10000.;
d = 1000.;

alpha_k = 100;
alpha_d = 100;

c  = get_c(q);
p = -get_ground_penetration(c);
J_c = get_J_c(q);
dc = -J_c* dq;
dp = get_ground_penetration(dc(2));

lambda = -k * (exp(alpha_k * p)-1) - d * sigmoid(alpha_d * p) * dp;

end

function s = sigmoid(x)
s = 1. / (1. + exp(-x));
end