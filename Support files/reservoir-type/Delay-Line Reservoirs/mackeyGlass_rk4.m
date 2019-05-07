%% mackeyGlass_rk4: Numerical solver of MG using 4-th order Runge Kutta
% This function computes the numerical solution of the Mackey-Glass delayed 
% differential equation using the 4-th order Runge-Kutta method

function x_t_plus_deltat = mackeyGlass_rk4(x_t, x_t_minus_tau, deltat, eta, T, p)

    k1 = deltat*mackeyGlass(x_t,          x_t_minus_tau, eta, T, p);
    k2 = deltat*mackeyGlass(x_t+0.5*k1,   x_t_minus_tau, eta, T, p);
    k3 = deltat*mackeyGlass(x_t+0.5*k2,   x_t_minus_tau, eta, T, p);
    k4 = deltat*mackeyGlass(x_t+k3,       x_t_minus_tau, eta, T, p);
    
    x_t_plus_deltat = (x_t + k1/6 + k2/3 + k3/3 + k4/6);

end