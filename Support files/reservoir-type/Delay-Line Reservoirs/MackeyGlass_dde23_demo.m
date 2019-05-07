%
% MackeyGlass_dde23_demo.m
%
% MATLAB demo script for the vector field: MackeyGlass
%
% This file was generated by the program VFGEN (Version:2.1.0)
% Generated on  3-Jun-2007 at 12:24
%
%
function MackeyGlass_dde23_demo(stoptime)
    a = 0.2;
    b = 0.1;
    tau = 17.0;
    p = zeros(3,1);
    p(1) = a;
    p(2) = b;
    p(3) = tau;
    p(4) = 10;
    
    % add input
    
    % timing parameters
    ts = 500; % start time
    T = 1000; % maximum elapsed time
    thold = 50;
    period = (T - ts) / thold;
    
    % set values
    x = zeros(1,1);
    x(1) = 0.5;

    lags = [tau];

    x0 = [0.5];

    opts = ddeset('reltol',1e-3,'abstol',1e-6,'InitialY',x0);
    sol = dde23(@(t,y,Z) MackeyGlass_dde23(t,y,Z,p),lags,@(t) MackeyGlass_history(t,p),[0 stoptime],opts);

    num_plot_samples = 500;
    tint = 1:num_plot_samples;%linspace(0,stoptime,num_plot_samples);
    xint = deval(sol,tint);
    
    clf
    plot(tint,xint,'linewidth',2);
    grid on
    xlabel('t');
    legend('x','Location','Best')
    
end

function x = MackeyGlass_history(t,p)
%     a          = p(1);
%     b          = p(2);
%     tau        = p(3);
    x(1) = 0.5+(0.02)*t;
end