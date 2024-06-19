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

% %% Centralized solution
[traj_central, xf_central] = central_sol(Planes, 0);

%% Dual decomposed solution
alpha = 0.5;
[~,xf_dual] = general_sol(Planes, opt_sim("dual", "constant", alpha, 0, 0, 0, "iteration", 3000, 1e-6), xf_central, 1);
logerr_plot(xf_central,xf_dual)

%% Investigate the effect of step size and step size update sequence
alpha = [0.2, 0.5];
xf_ac = cell(1, length(alpha));
for i = 1:length(alpha)
    xf_ac{i} = general_sol(Planes, opt_sim("dual", "constant", alpha(i), 0, 0, 0, "iteration", 3000, 1e-6), xf_central, 1);
end
a1e_1_var = general_sol(Planes, opt_sim("dual", "variable", 5e-1, gamma, rho, phi, "iteration", 3000, 1e-6), xf_central, 1);

figure(34), clf;
hold on
plot(xf_ac{1});
plot(a1e_1_var);
yscale('log')

%% Nesterov accelerated method of subgradient updates
a1e_1_con = general_sol(Planes, opt_sim("dual", "constant", 4e-1, 0.8, rho, phi, "iteration", 3000, 1e-6), xf_central, 1);
a1e_1_var = general_sol(Planes, opt_sim("dual", "variable", 9e-1, 0.8, rho, phi, "iteration", 3000, 1e-6), xf_central, 1);
a1e_1_nes = general_sol(Planes, opt_sim("dual", "nesterov", 3e-1, 0.8, rho, phi, "iteration", 3000, 1e-6), xf_central, 1);

figure(34), clf;
hold on
plot(a1e_1_con);
plot(a1e_1_var);
plot(a1e_1_nes);
yscale('log')

%% Consensus based approach
phi = [1, 4, 10, 100];
alpha = [0.4, 0.4, 0.4, 0.4];
xf_phi = cell(1, length(phi));
for i = 1:length(phi)
    xf_phi{i} = general_sol(Planes, opt_sim("consensus", "constant", 0.9, 0, 0, phi(i), "iteration", 5000, 1e-6), xf_central, 0);
end

figure(34), clf;
tiledlayout(2,2);
hold on
for i = 1:length(phi)
    nexttile
    plot(xf_phi{i}')
    yscale('log')
    ylim([1e-4, 1e1])
    xlim([0,5000])
end
% legend

%% Consensus ADMM
rho = [1.5, 2, 5, 10];
xf_rho = cell(1, length(rho));
for i = 1:length(rho)
    xf_rho{i} = general_sol(Planes, opt_sim("ADMM", "constant", 0, 0, rho(i), 0, "iteration", 1000, 1e-6), xf_central, 1);
end

figure(34), clf;
hold on
for i = 1:length(rho)
    plot(xf_rho{i}')
end
yscale('log')
legend
