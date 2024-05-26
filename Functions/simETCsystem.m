function [tTotal, xTotal, events] = simETCsystem(tspan, x0, sigma, A, B, K, P, Q)
    % Initialize system
    tTotal = tspan(1); %Starting time
    xi_0 = x0;
    eps_0 = [0;0];
    xTotal = [xi_0;eps_0]'; %initialize with zero error
    events = 0;

    Efun = @(t,x) eventFunction(t, x, sigma, B, K, P, Q);

    % Run simulation
    while tTotal(end) < tspan(2)
        options = odeset('Events', @(t,x) Efun(t,x), 'RelTol', 1e-5, 'AbsTol', 1e-7);
        odeFun = @(t, xi) systemDynamics(t, xi, A, B, K);
        tRemainder = [tTotal(end), tspan(2)]; % Remainder of Tspan
        xReset = [xTotal(end, 1:2),0,0]'; % Reset the error back to zero
        [t, x, te, ~, ~] = ode45(odeFun, tRemainder, xReset, options);
        tTotal = [tTotal; t(2:end)];
        xTotal = [xTotal; x(2:end, :)];
        if ~isempty(te)
            if size(te) ~= 1
                disp(size(te))
            end
            events = events + 1; % Event occurred
        else
            break; % Simulation is finished
        end
    end
end