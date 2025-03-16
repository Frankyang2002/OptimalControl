function x_next = rk4_step(f, x, dt)
    % f: a function handle that takes a single argument, x
    % x: current state
    % dt: time step

    k1 = f(x);
    k2 = f(x + dt/2 * k1);
    k3 = f(x + dt/2 * k2);
    k4 = f(x + dt   * k3);

    x_next = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end
