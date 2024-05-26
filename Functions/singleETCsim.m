function singleETCsim(sigma, x0, T_end, dist)

    a = 5;
    b = 9;
    c = 8;
    
    A = [0.3+a-b, 0.5-c;
         0, 1];
    B = [0;1];
    
    p = [-1-2j, -1+2j];
    K = place(A,B,p);
    
    Q = eye(2);
    cvx_begin sdp quiet
        variable P(2,2) semidefinite
        subject to
            P == P.'
            (A-B*K)'*P + P*(A-B*K) == -Q
    cvx_end

    % Simulation parameters
    tspan = [0 T_end];
    
    % Initialize zero matrices
    simETC = @(tspan, x0, sigma) simETCsystem(tspan, x0, sigma, A, B, K, P, Q, dist);
    
    tic
    [tTotal, xTotal, ~] = simETC(tspan, x0, sigma);
    figure(dist),clf;
    tiledlayout(2,1);
    nexttile
    plot(tTotal, xTotal(:,1:2));
    grid on
    nexttile
    plot(tTotal, xTotal(:,3:4));
    grid on
    toc
end