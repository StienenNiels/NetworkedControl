function singleETCsim(sigma, x0, T_end, dist, eps)

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
    simETC = @(tspan, x0, sigma) simETCsystem_eps(tspan, x0, sigma, A, B, K, P, Q, dist, eps);

    % Figure label generation
    if dist == 0
        lgd1 = "$d(t) = [0,0]^T$";
    elseif dist == 1
        lgd1 = "$d(t) = 0.1[sin(t), cos(t)]^T$"; 
    elseif dist == 2
        lgd1 = "$d(t) = [1,1]^Trnd()$";
    else
        lgd1 = "$d(t) = 0.1[|sin(t)|, -|cos(t)|]^T$"; 
    end

    if eps == 0
        title_str = "$\sigma = " + num2str(sigma) + "$ and " + lgd1;
    else
        title_str = "$\sigma = " + num2str(sigma) + "$ and " + lgd1 + " and $\epsilon = " + num2str(eps) + "$";
    end

    tic
    [tTotal, xTotal, ~] = simETC(tspan, x0, sigma);
    fig = figure(99);
    clf;
    tiledlayout(2,1);
    nexttile
    plot(tTotal, xTotal(:,1:2));
    xlabel("time [seconds]")
    ylabel("magnitude")
    legend("$\xi_1$", "$\xi_2$");
    title(title_str);
    grid on
    nexttile
    plot(tTotal, xTotal(:,3:4));
    xlabel("time [seconds]")
    ylabel("magnitude")
    legend("$\epsilon_1$", "$\epsilon_2$");
    grid on
    toc
    set(fig, 'Position', [1000, 100, 500, 400]);
    set(gcf, "Theme", "light");
    filename = "Images\2_491_S" + num2str(round(100*sigma)) + "_D" + num2str(dist) + "_E" + num2str(strlength(num2str(eps))) + ".pdf";
    exportgraphics(gcf,filename)
    close(gcf);
end