function bulk_hs_sim(Ninit,init_max, T_end, hs, dist)
    a = 5;
    b = 9;
    c = 8;
    
    A = [0.3+a-b, 0.5-c;
         0, 1];
    B = [0;1];
    
    p = [-1-2j, -1+2j];
    K = place(A,B,p);

    % Simulation parameters
    x0_set = -init_max + init_max * rand(2, Ninit);
    tspan = [0 T_end];
    
    % Initialize zero matrices
    stable = zeros(1, size(x0_set, 2));
    sim_hs = @(tspan, x0) sim_hs_system(tspan, x0, A, B, K, hs, dist);
    
    tic
    parfor init_cond_ind = 1:size(x0_set, 2)
        x0 = x0_set(:,init_cond_ind);
        [~, xTotal, ~] = sim_hs(tspan, x0);
        if round(xTotal(end,1:2)) == 0 % Check if solution was stable
            stable(init_cond_ind) = 1;
        else
            disp([init_cond_ind, round(xTotal(end,1:2))])
        end
    end
    toc
    
    p_stable = 100*sum(sum(stable))/(size(x0_set, 2));
    fprintf("Percentage of results that are stable: %.2f%%\n", p_stable)
end