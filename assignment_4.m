clear
clc

addpath("Functions_A4\")
addpath("Images\")

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% state vector is (x,y,xdot,ydot)
aircraft; % Load the given parameters
Planes = planes_gen(); % Define structures

% Create structures for each plane to make workspace more structured
for i = 1:4
    Planes(i).plane = i;
    Planes(i).Tf = Tfinal;
    Planes(i).umax = umax;
    Planes(i).A = eval(sprintf('A%d;', i));
    Planes(i).B = eval(sprintf('B%d;', i));
    Planes(i).x0 = eval(sprintf('x0%d;', i));
    Planes(i) = predmodgen(Planes(i));
    Planes(i) = optgen(Planes(i));
end
dim = Planes(1).dim;
clearvars -except Planes Tfinal umax dim


%% Centralized solution
[traj_central, xf_central] = central_sol(Planes, 0);

%% Dual decomposed solution
alpha = 0.1;
[traj_dual, xf_dual] = dual_sol(Planes, alpha, 0, 0); % 0 for constant alpha
%%
% Still need to make the figures of the trajectories
%
logerr_plot(xf_central,xf_dual)

%% Investigate the effect of step size and step size update sequence
[~, a1e_1] = dual_sol(Planes, 1e-1, 0, 1); % 0 for constant alpha
avg_a1e_1 = err_norm(xf_central,a1e_1,1);
[~, a1e_1_var] = dual_sol(Planes, 5e-1, 1, 1); % 1 for variable alpha
avg_a1e_1_var = err_norm(xf_central,a1e_1_var,1);

figure(34), clf;
hold on
plot(avg_a1e_1);
plot(avg_a1e_1_var);
yscale('log')
% [~, a1e_3] = dual_sol(Planes, 1e-3, 0, 0); % 0 for constant alpha
% [~, a1e_4] = dual_sol(Planes, 1e-4, 0, 0); % 0 for constant alpha


%% Nesterov accelerated method of subgradient updates
[~, a1e_1] = dual_sol(Planes, 4e-1, 0, 0); % 0 for constant alpha
avg_a1e_1 = err_norm(xf_central,a1e_1,1);
[~, a1e_1_var] = dual_sol(Planes, 9e-1, 1, 0); % 1 for variable alpha
avg_a1e_1_var = err_norm(xf_central,a1e_1_var,1);
[~, a1e_1_nes] = dual_sol(Planes, 1e-1, 2, 0); % 2 for nesterov's
avg_a1e_1_nes = err_norm(xf_central,a1e_1_nes,1);

figure(34), clf;
hold on
plot(avg_a1e_1);
plot(avg_a1e_1_var);
plot(avg_a1e_1_nes);
yscale('log')

% Just need to make plots of the results

%% Consensus based approach
[~, a1e_1] = consensus_sol(Planes, 4e-1, 0, 0);
avg_a1e_1 = err_norm(xf_central,a1e_1,1);
[~, a1e_2] = consensus_sol(Planes, 4e-1, 1, 0);
avg_a1e_2 = err_norm(xf_central,a1e_2,1);
[~, a1e_3] = consensus_sol(Planes, 4e-1, 2, 0);
avg_a1e_3 = err_norm(xf_central,a1e_3,1);
[~, a1e_4] = consensus_sol(Planes, 4e-1, 5, 0);
avg_a1e_4 = err_norm(xf_central,a1e_4,1);
[~, a1e_5] = consensus_sol(Planes, 4e-1, 10, 0);
avg_a1e_5 = err_norm(xf_central,a1e_5,1);
[~, a1e_6] = consensus_sol(Planes, 4e-1, 20, 0);
avg_a1e_6 = err_norm(xf_central,a1e_6,1);

figure(34), clf;
hold on
plot(avg_a1e_1);
plot(avg_a1e_2);
plot(avg_a1e_3);
plot(avg_a1e_4);
plot(avg_a1e_5);
plot(avg_a1e_6);
yscale('log')