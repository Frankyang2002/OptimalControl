function p = get_ground_penetration(x)

%p= 0.5*tanh(-10000.*x).*x +0.5*x;
p= 0.5*sign(-x) *x + 0.5*x;

end