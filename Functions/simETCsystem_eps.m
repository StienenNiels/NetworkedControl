function [tTotal, xTotal, events] = simETCsystem_eps(tspan, x0, sigma, A, B, K, P, Q, dist)
    % Initialize system
    Nguess = round(tspan(2)/0.01);
    tTotal = zeros(1, Nguess);
    xTotal = zeros(Nguess, 4);
    tTotal(1) = tspan(1); %Starting time
    xi_0 = x0;
    eps_0 = [0;0];
    xTotal(1,:) = [xi_0;eps_0]'; %initialize with zero error
    events = 0;

    ind = 1;

    Efun = @(t,x) eventFunction_eps(t, x, sigma, B, K, P, Q);

    % Run simulation
    while tTotal(ind) < tspan(2)
        options = odeset('Events', @(t,x) Efun(t,x), 'RelTol', 1e-5, 'AbsTol', 1e-7);
        odeFun = @(t, xi) systemDynamics(t, xi, A, B, K, dist);
        tRemainder = [tTotal(ind), tspan(2)]; % Remainder of Tspan
        xReset = [xTotal(ind, 1:2),0,0]'; % Reset the error back to zero
        [t, x, te, ~, ~] = ode45(odeFun, tRemainder, xReset, options);
        
        Nnew = length(t) - 1;
        if ind + Nnew > length(tTotal)
            % Double the size of the arrays
            tTotal = [tTotal, zeros(1, length(tTotal))];
            xTotal = [xTotal; zeros(size(xTotal))];
        end
        tTotal(ind + 1 : ind + Nnew) = t(2:end);
        xTotal(ind + 1 : ind + Nnew, :) = x(2:end, :);
        ind = ind + Nnew;

        if ~isempty(te)
            if size(te) ~= 1
                disp(size(te))
            end
            events = events + 1; % Event occurred
        else
            break; % Simulation is finished
        end
    end
    tTotal = tTotal(1:ind);
    xTotal = xTotal(1:ind, :);
end