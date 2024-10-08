clear
clc

addpath("Functions\")
addpath("Images\")

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
%%

% Number of sigma values to test:
Nsigma = 10;
sigma_min = 0.01;
sigma_max = 0.99;
% Number of initial conditions to test
Ninit = 2;
init_max = 5;

T_end = 10;
plotastate = 0;

%%
Ninit = 5;
[sig_val, total_events] = bulkETCsim(10,0.01,0.95, Ninit,2, 10, 0);
hs = T_end./(total_events/Ninit);
[~,hs_ind] = max(hs);

% figure(411), clf;
% colororder({'k','#0072BD'})
% 
% yyaxis left;
% hold on;
% plot(sig_val, hs, '-o', 'LineWidth', 1.5);
% plot(sig_val(hs_ind), hs(hs_ind),'.', "MarkerSize",25, "Color","#77AC30");
% ylim([0, 0.12])
% yll = ylabel('$h_{avg}$');
% fontsize(yll,15,"points");
% 
% yyaxis right;
% plot(sig_val, total_events/Ninit, '-x', 'LineWidth', 1.3, "MarkerSize",10);
% ylim([0, ceil(max(total_events/Ninit)/6000)*6000])
% ylr = ylabel('$count(s_k)$');
% fontsize(ylr,15,"points");
% 
% xl = xlabel('$\sigma$');
% fontsize(xl,15,"points");
% xlim([0,1])
% lgd = legend('$h_{avg}$', ['$h^*_{avg} = ', num2str(hs(hs_ind), '%.3f'),'$'], '$avg(N_{comms})$', ...
%      "Location","northeast");
% fontsize(lgd,11,"points");
% grid on;
% set(gcf, "Theme", "light"); % Uncomment for report plots

%%
Ninit = 5;
[sig_val, total_events] = bulkETCsim(11,0.01,0.3, Ninit,2, 10, 0);
hs = T_end./(total_events/Ninit);
[~,hs_ind] = max(hs);

% figure(412), clf;
% colororder({'k','#0072BD'})
% 
% yyaxis left;
% hold on;
% plot(sig_val, hs, '-o', 'LineWidth', 1.5);
% plot(sig_val(hs_ind), hs(hs_ind),'.', "MarkerSize",25, "Color","#77AC30");
% ylim([0, 0.12])
% yll = ylabel('$h_{avg}$');
% fontsize(yll,15,"points");
% 
% yyaxis right;
% plot(sig_val, total_events/Ninit, '-x', 'LineWidth', 1.3, "MarkerSize",10);
% ylim([0, ceil(max(total_events/Ninit)/500)*500])
% ylr = ylabel('$count(s_k)$');
% fontsize(ylr,15,"points");
% 
% xl = xlabel('$\sigma$');
% fontsize(xl,15,"points");
% xlim([0,0.3])
% lgd = legend('$h_{avg}$', ['$h^*_{avg} = ', num2str(hs(hs_ind), '%.3f'),'$'], '$avg(N_{comms})$', ...
%      "Location","northeast");
% fontsize(lgd,11,"points");
% grid on;
% set(gcf, "Theme", "light"); % Uncomment for report plots

%% Question 4.3
sigma = 0.039;
hs_avg = 0.104;
Ninit = 10;

bulk_hs_sim(Ninit,2, 10, hs_avg, 0)


a = 5;
b = 9;
c = 8;

A = [0.3+a-b, 0.5-c;
     0, 1];
B = [0;1];

p = [-1-2j, -1+2j];
K = place(A,B,p);
A_cl = c2d_zoh(A, B, K, 0, hs_avg, 0);
stable_43 = sr(A_cl)

%% Question 4.4
clc
Ninit = 5;
init_max = 1;
x0_set = -init_max + init_max * rand(2, Ninit);
T_end = 20;
[sig_val, total_events0] = bulkETCsim_dist(10,0.01,0.95, x0_set, T_end, 0);
[~, total_events1]       = bulkETCsim_dist(10,0.01,0.95, x0_set, T_end, 1);
[~, total_events2]       = bulkETCsim_dist(10,0.01,0.95, x0_set, T_end, 2);
[~, total_events3]       = bulkETCsim_dist(10,0.01,0.95, x0_set, T_end, 3);
hs0 = T_end./(total_events0/Ninit);
hs1 = T_end./(total_events1/Ninit);
hs2 = T_end./(total_events2/Ninit);
hs3 = T_end./(total_events3/Ninit);
% [~,hs_ind] = max(hs1);

figure(441), clf;
hold on;
plot(sig_val, hs0, '-o', 'LineWidth', 1.5);
plot(sig_val, hs1, '-o', 'LineWidth', 1.5);
plot(sig_val, hs2, '-o', 'LineWidth', 1.5);
plot(sig_val, hs3, '-o', 'LineWidth', 1.5);
% plot(sig_val(hs_ind), hs(hs_ind),'.', "MarkerSize",25, "Color","#77AC30");
% ylim([0, 0.05])
yll = ylabel('$h_{avg}$');
fontsize(yll,15,"points");
xl = xlabel('$\sigma$');
fontsize(xl,15,"points");
xlim([0,1])
lgd = legend('$d_0=[0,0]^T$','$d_1=0.1[sin(t), cos(t)]^T$', '$d_2=[1,1]^Trnd()$', '$d_3=0.1[|sin(t)|, -|cos(t)|]^T$', ...
     "Location","northeast");
fontsize(lgd,11,"points");
grid on;
set(gcf, "Theme", "light"); % Uncomment for report plots



%% Question 4.5
clc
Ninit = 10;
init_max = 1;
x0_set = -[init_max;init_max] + init_max * rand(2, Ninit);
T_end = 20;
[sig_val, total_events0] = bulkETCsim_dist(10,0.01,0.95, x0_set, T_end, 0);
[~, total_events1]       = bulkETCsim_eps(10,0.01,0.95, x0_set, T_end, 0, 1e-8);
[~, total_events2]       = bulkETCsim_eps(10,0.01,0.95, x0_set, T_end, 0, 1e-6);
[~, total_events3]       = bulkETCsim_eps(10,0.01,0.95, x0_set, T_end, 0, 1e-4);
hs0 = T_end./(total_events0/Ninit);
hs1 = T_end./(total_events1/Ninit);
hs2 = T_end./(total_events2/Ninit);
hs3 = T_end./(total_events3/Ninit);
% [~,hs_ind] = max(hs1);

figure(451), clf;
hold on;
plot(sig_val, hs0, '-o', 'LineWidth', 1.5);
plot(sig_val, hs1, '-o', 'LineWidth', 1.5);
plot(sig_val, hs2, '-o', 'LineWidth', 1.5);
plot(sig_val, hs3, '-o', 'LineWidth', 1.5);
% plot(sig_val(hs_ind), hs(hs_ind),'.', "MarkerSize",25, "Color","#77AC30");
% ylim([0, 0.05])
yll = ylabel('$h_{avg}$');
fontsize(yll,15,"points");
xl = xlabel('$\sigma$');
fontsize(xl,15,"points");
xlim([0,1])
lgd = legend('$\phi(\xi(t), \xi(s_k)) \leq 0$','$\phi(\xi(t), \xi(s_k)) \leq 10^{-8}$','$\phi(\xi(t), \xi(s_k)) \leq 10^{-6}$','$\phi(\xi(t), \xi(s_k)) \leq 10^{-4}$', ...
     "Location","northeast");
fontsize(lgd,11,"points");
grid on;
set(gcf, "Theme", "light"); % Uncomment for report plots


%% Generate some state trajectories
clc

% State trajectories for sigma variations
singleETCsim(0.01, [1;0], 10, 0, 0)
singleETCsim(0.4, [1;0], 10, 0, 0)
singleETCsim(0.8, [1;0], 10, 0, 0)

%%
% State trajectories for ic variations
singleETCsim(0.04, [10;3], 10, 0, 0)
singleETCsim(0.04, [-5;7], 10, 0, 0)
singleETCsim(0.04, [1;-1], 10, 0, 0)

%%
% State trajectories for disturbance variations
singleETCsim(0.04, [1;0], 30, 0, 0)
singleETCsim(0.04, [1;0], 30, 1, 0)
singleETCsim(0.04, [1;0], 30, 2, 0)
singleETCsim(0.04, [1;0], 30, 3, 0)

%%
% State trajectories for epsilon variations
singleETCsim(0.04, [1;0], 10, 0, 0)
singleETCsim(0.04, [1;0], 10, 0, 1e-8)
singleETCsim(0.04, [1;0], 10, 0, 1e-6)
singleETCsim(0.04, [1;0], 10, 0, 1e-4)
singleETCsim(0.04, [1;0], 10, 0, 1e-2)
singleETCsim(0.04, [1;0], 10, 0, 0.9)


