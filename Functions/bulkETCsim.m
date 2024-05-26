function [sig_val, total_events] = bulkETCsim(Nsigma,sigma_min,sigma_max, Ninit,init_max, T_end, dist)

    student_id = 5595738;
    a = 5;
    b = 9;
    c = 8;
    
    A = [0.3+a-b, 0.5-c;
         0, 1];
    B = [0;1];
    
    p = [-1-2j, -1+2j];
    K = place(A,B,p);
    
    lambdas = eig(A-B*K);
    
    Q = eye(2);
    cvx_begin sdp quiet
        variable P(2,2) semidefinite
        subject to
            P == P.'
            (A-B*K)'*P + P*(A-B*K) == -Q
    cvx_end

    % Simulation parameters
    sig_val = linspace(sigma_min,sigma_max, Nsigma);
    x0_set = -init_max + init_max * rand(2, Ninit);
    tspan = [0 T_end];
    
    % Initialize zero matrices
    total_events = zeros(length(sig_val),1);
    stable = zeros(length(sig_val), size(x0_set, 2));
    simETC = @(tspan, x0, sigma) simETCsystem(tspan, x0, sigma, A, B, K, P, Q, dist);
    
    tic
    for ind = 1:length(sig_val)
        sigma = sig_val(ind);
        event_tot = 0;
        parfor init_cond_ind = 1:size(x0_set, 2)
            x0 = x0_set(:,init_cond_ind);
            [~, xTotal, events] = simETC(tspan, x0, sigma);
            event_tot = event_tot + events;
            if dist == 0
                if round(xTotal(end,1:2)) == 0 % Check if solution was stable
                    stable(ind, init_cond_ind) = 1;
                else
                    disp([sigma, init_cond_ind, round(xTotal(end,1:2))])
                end
            else
                if round(xTotal(end,1:2)/10)*10 == 0 % Check if solution was stable
                    stable(ind, init_cond_ind) = 1;
                else
                    disp([sigma, init_cond_ind, round(xTotal(end,1:2))])
                end
            end
        end
        % if plotastate
        %     figure(1),clf;
        %     plot(tTotal, xTotal);
        %     grid on
        % end
        total_events(ind) = event_tot;
    end
    toc
    
    p_stable = 100*sum(sum(stable))/(length(sig_val)*size(x0_set, 2));
    fprintf("Percentage of results that are stable: %.2f%%\n", p_stable)
end