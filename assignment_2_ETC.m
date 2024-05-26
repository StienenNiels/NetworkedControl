clear
clc

addpath("Functions\")

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%% Initialization of matrices A and B
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
%%

% Number of sigma values to test:
Nsigma = 10;
% Number of initial conditions to test
Ninit = 10;

% Simulation parameters
sig_val = linspace(0.01,0.9, Nsigma);
x0_set = -1 + 2 * rand(2, Ninit);
tspan = [0 10];

% Initialize zero matrices
total_events = zeros(length(sig_val),1);
event_matrix = zeros(length(sig_val), size(x0_set, 2));
stable = zeros(length(sig_val), size(x0_set, 2));
simETC = @(tspan, x0, sigma) simETCsystem(tspan, x0, sigma, A, B, K, P, Q);

tic
for ind = 1:length(sig_val)
    sigma = sig_val(ind);
    event_tot = 0;
    for init_cond_ind = 1:size(x0_set, 2)
        x0 = x0_set(:,init_cond_ind);
        [tTotal, xTotal, events] = simETC(tspan, x0, sigma);
        event_matrix(ind,init_cond_ind) = events;
        event_tot = event_tot + events;
        if round(xTotal(end,1:2),1) == 0 % Check if solution was stable
            stable(ind, init_cond_ind) = 1;
        else
            disp([sigma, init_cond_ind, round(xTotal(end,1:2),2)])
        end
    end
    total_events(ind) = event_tot;
end
toc

p_stable = 100*sum(sum(stable))/(length(sig_val)*size(x0_set, 2));
fprintf("Percentage of results that are stable: %.2f%%\n", p_stable)

% Plotting the results
figure(4),clf;
hold on;
plot(sig_val, total_events/Ninit, '-o');
xlabel('Sigma Values');
ylabel('Average Number of Communications');
title('Communications vs. Sigma');
xlim([0,1])
ylim([0, ceil(max(total_events/Ninit)/1000)*1000])
grid on;

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
            events = events + 1; % Event occurred
        else
            break; % Simulation is finished
        end
    end
end

function [value, isterminal, direction] = eventFunction(~, xi_eps, sigma, B, K, P, Q)
    value = xi_eps'*[(1-sigma)*Q, P*B*K; K'*B'*P, zeros(2)]*xi_eps; % Slide 49/57
    isterminal = 1; % Stop the integration
    direction = 0; % Detect all zero crossings
end

function xi_eps_dot = systemDynamics(~, xi_eps, A, B, K) % Slide 48/57
    xi = xi_eps(1:2);
    eps = xi_eps(3:4);

    xi_dot = (A-B*K)*xi -B*K*eps;
    eps_dot = -(A-B*K)*xi +B*K*eps;
    xi_eps_dot = [xi_dot; eps_dot];    
end