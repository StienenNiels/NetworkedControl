function xi_eps_dot = systemDynamics(t, xi_eps, A, B, K, dist) % Slide 48/57
    xi = xi_eps(1:2);
    eps = xi_eps(3:4);

    if dist == 0
        xi_dot = (A-B*K)*xi -B*K*eps;
    elseif dist == 1
        xi_dot = (A-B*K)*xi -B*K*eps + 1*[1;1]*sin(t);
    elseif dist == 2
        xi_dot = (A-B*K)*xi -B*K*eps + 2*[1;1]*sin(t);
    elseif dist == 3
        xi_dot = (A-B*K)*xi -B*K*eps + 3*[1;1]*sin(t);
    end
    eps_dot = -(A-B*K)*xi +B*K*eps;
    xi_eps_dot = [xi_dot; eps_dot];    
end