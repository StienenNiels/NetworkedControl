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
clearvars -except Planes

% Centralized solution
[traj_central, xf_central] = central_sol(Planes);

%% Dual decomposed solution
alpha = 0.5;
[~,xf_dual,traj_dual] = general_sol(Planes, opt_sim("dual", "constant", alpha, 0, 0, 0, "iteration", 2000, 1e-6), xf_central, 1);
% 3 plots needed
logerr_plot(xf_central,xf_dual, "1a_error")
xyplane_plot(xf_central,traj_dual, "1a_XY")
trajectories_plot(xf_central,traj_dual, "1a_traj")
clearvars -except Planes xf_central

%% Investigate the effect of step size and step size update sequence
alpha = [0.1, 0.2, 0.4, 0.6, 0.9];
xf_ac = cell(1, length(alpha));
xf_av1 = cell(1, length(alpha));
xf_av2 = cell(1, length(alpha));
xf_av3 = cell(1, length(alpha));
for i = 1:length(alpha)
    xf_ac{i}  = general_sol(Planes, opt_sim("dual", "constant",  alpha(i), 0, 0, 0, "iteration", 2000, 1e-6), xf_central, 1);
    xf_av1{i} = general_sol(Planes, opt_sim("dual", "variable1", alpha(i), 0, 0, 0, "iteration", 2000, 1e-6), xf_central, 1);
    xf_av2{i} = general_sol(Planes, opt_sim("dual", "variable2", alpha(i), 0, 0, 0, "iteration", 2000, 1e-6), xf_central, 1);
    xf_av3{i} = general_sol(Planes, opt_sim("dual", "variable3", alpha(i), 0, 0, 0, "iteration", 2000, 1e-6), xf_central, 1);
end
plot_info.alpha = alpha;
logerr_cell_plot(xf_ac, plot_info, "constant", "1b_constant")
logerr_cell_plot(xf_av1, plot_info, "variable1", "1b_variable1")
logerr_cell_plot(xf_av2, plot_info, "variable2", "1b_variable2")
logerr_cell_plot(xf_av3, plot_info, "variable3", "1b_variable3")
clearvars -except Planes xf_central

%% Nesterov accelerated method of subgradient updates
alpha = [0.1, 0.1, 0.1, 0.9, 0.9, 0.9];
gamma = [0.5, 0.8, 0.9, 0.5, 0.8, 0.9];
xf_nes = cell(1, length(alpha));
for i = 1:length(alpha)
    xf_nes{i}  = general_sol(Planes, opt_sim("dual", "nesterov",  alpha(i), gamma(i), 0, 0, "iteration", 500, 1e-6), xf_central, 1);
end

plot_info.alpha = alpha;
plot_info.gamma = gamma;
logerr_cell_plot(xf_nes, plot_info, "nesterov", "1c_nesterov")
clearvars -except Planes xf_central

%% Comparison of all the above methods
alpha = [0.7, 0.9, 0.9, 0.9, 0.2];
gamma = [0.0, 0.0, 0.0, 0.0, 0.8];
xf_comp = cell(1, 5);
xf_comp{1} = general_sol(Planes, opt_sim("dual", "constant",  alpha(1), gamma(1), 0, 0, "iteration", 2000, 1e-6), xf_central, 1);
xf_comp{2} = general_sol(Planes, opt_sim("dual", "variable1", alpha(2), gamma(2), 0, 0, "iteration", 2000, 1e-6), xf_central, 1);
xf_comp{3} = general_sol(Planes, opt_sim("dual", "variable2", alpha(3), gamma(3), 0, 0, "iteration", 2000, 1e-6), xf_central, 1);
xf_comp{4} = general_sol(Planes, opt_sim("dual", "variable3", alpha(4), gamma(4), 0, 0, "iteration", 2000, 1e-6), xf_central, 1);
xf_comp{5} = general_sol(Planes, opt_sim("dual", "nesterov",  alpha(5), gamma(5), 0, 0, "iteration", 2000, 1e-6), xf_central, 1);

plot_info.alpha = alpha;
plot_info.gamma = gamma;
logerr_cell_plot(xf_comp, plot_info, "dual_comp", "1c_comparison")
clearvars -except Planes xf_central

%% Consensus based approach
phi = [1, 5, 20, 50, 10, 50];
alpha = [0.4, 0.4, 0.4, 0.4, 0.7, 0.7];
xf_phi = cell(1, length(phi));
for i = 1:length(phi)
    [xf_phi{i},xf] = general_sol(Planes, opt_sim("consensus", "constant", alpha(i), 0, 0, phi(i), "iteration", 1500, 1e-6), xf_central, 1);
end

plot_info.alpha = alpha;
plot_info.phi = phi;
logerr_cell_plot(xf_phi, plot_info, "consensus", "1d_consensus")
clearvars -except Planes xf_central

%% Consensus ADMM
rho = 5;
[~,xf_ADMM,traj_ADMM] = general_sol(Planes, opt_sim("ADMM", "constant", 0, 0, rho, 0, "iteration", 200, 1e-6), xf_central, 1);
logerr_plot(xf_central,xf_ADMM, "2a_error")
xyplane_plot(xf_central,traj_ADMM, "2a_XY")
trajectories_plot(xf_central,traj_ADMM, "2a_traj")
clearvars -except Planes xf_central

%% Consensus ADMM varying rho
rho = [1, 1.5, 2, 5, 6.5, 8.5, 10.5, 20];
xf_rho = cell(1, length(rho));
for i = 1:length(rho)
    xf_rho{i} = general_sol(Planes, opt_sim("ADMM", "constant", 0, 0, rho(i), 0, "iteration", 500, 1e-6), xf_central, 1);
end

plot_info.rho = rho;
logerr_cell_plot(xf_rho, plot_info, "ADMM", "2b_rho")
clearvars -except Planes xf_central

%% Consensus ADMM optimal rho
% Needs to be balanced between fast convergence and lower steady state error
rho = 0.1:0.1:50;
xf_rho = cell(1, length(rho));
for i = 1:length(rho)
    xf_rho{i} = general_sol(Planes, opt_sim("ADMM", "constant", 0, 0, rho(i), 0, "convergence", 220, 1e-5), xf_central, 1);
    if mod(i, 20) == 0
        disp(i);
    end
end

plot_info.rho = rho;
rho_opt_plot(xf_rho, plot_info, "rho", "2b_rho_optimal")
clearvars -except Planes xf_central