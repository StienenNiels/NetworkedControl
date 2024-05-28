function [value, isterminal, direction] = eventFunction_eps(~, xi_eps, sigma, B, K, P, Q)
    eps = 0.0000001;
    value = xi_eps'*[(1-sigma)*Q, P*B*K; K'*B'*P, zeros(2)]*xi_eps + eps; % Slide 49/57
    isterminal = value <= 0; % Stop the integration
    direction = 0; % Detect all zero crossings
end