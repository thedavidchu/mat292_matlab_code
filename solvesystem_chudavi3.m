%2D IEM
%solvesystem_chudavi3

function [t_out,x_out] = solvesystem_chudavi3(f,g,time_interval,x_0,h)
    %Time interval: [t0,tN]
    t_0 = time_interval(1);
    t_N = time_interval(2);

    %Step size: h
    N = ceil((t_N - t_0)/h); %Check endpoints
    t_out = linspace(t_0,t_N,N);

    %Initial condition: 
    x_a = zeros(N,1);
    x_b = zeros(N,1);

    x_a(1) = x_0(1);
    x_b(1) = x_0(2);

    for i = 1:N-1
        %Use vectors instead of all this junk!
        %Equation for Euler's method y = y_prev + slope * (t - t_prev);
        x_a(i+1) = x_a(i) + h/2*( f(t_out(i),x_a(i),x_b(i)) + f(t_out(i)+h, x_a(i)+f(t_out(i),x_a(i),x_b(i))*h, x_b(i) + g(t_out(i),x_a(i),x_b(i))*h));
        x_b(i+1) = x_b(i) + h/2*( g(t_out(i),x_a(i),x_b(i)) + g(t_out(i)+h, x_a(i)+f(t_out(i),x_a(i),x_b(i))*h, x_b(i) + g(t_out(i),x_a(i),x_b(i))*h));
    end
    
    x_out = transpose([x_a, x_b]);
end