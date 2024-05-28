function [value, isterminal, direction] = eventFunction(~, xi_eps, sigma, B, K, P, Q)
    value = xi_eps'*[(1-sigma)*Q, P*B*K; K'*B'*P, zeros(2)]*xi_eps; % Slide 49/57
    isterminal = value <= 0; % Stop the integration
    direction = 0; % Detect all zero crossings
end