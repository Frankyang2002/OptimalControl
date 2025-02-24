function x_next = rk4_step(f, x, dt)
    import casadi.*;
    k1 = f(x);
    k2 = f(x + dt/2*k1);
    k3 = f(x + dt/2*k2);
    k4 = f(x + dt*k3);
    x_next = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);

end