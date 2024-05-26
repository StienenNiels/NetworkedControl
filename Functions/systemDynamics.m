function xi_eps_dot = systemDynamics(~, xi_eps, A, B, K) % Slide 48/57
    xi = xi_eps(1:2);
    eps = xi_eps(3:4);

    xi_dot = (A-B*K)*xi -B*K*eps;
    eps_dot = -(A-B*K)*xi +B*K*eps;
    xi_eps_dot = [xi_dot; eps_dot];    
end