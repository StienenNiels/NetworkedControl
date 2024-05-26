clear
clc

addpath("Functions\")

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

%% Question 4.4
% clc
% Ninit = 5;
% [sig_val, total_events1] = bulkETCsim(10,0.01,0.95, Ninit,1, 30, 1);
% [~, total_events2]       = bulkETCsim(10,0.01,0.95, Ninit,1, 30, 2);
% [~, total_events3]       = bulkETCsim(10,0.01,0.95, Ninit,1, 30, 3);
% hs1 = T_end./(total_events1/Ninit);
% hs2 = T_end./(total_events2/Ninit);
% hs3 = T_end./(total_events3/Ninit);
% [~,hs_ind] = max(hs1);
% 
% figure(441), clf;
% hold on;
% plot(sig_val, hs1, '-o', 'LineWidth', 1.5);
% plot(sig_val, hs2, '-o', 'LineWidth', 1.5);
% plot(sig_val, hs3, '-o', 'LineWidth', 1.5);
% % plot(sig_val(hs_ind), hs(hs_ind),'.', "MarkerSize",25, "Color","#77AC30");
% ylim([0, 0.05])
% yll = ylabel('$h_{avg}$');
% fontsize(yll,15,"points");
% xl = xlabel('$\sigma$');
% fontsize(xl,15,"points");
% xlim([0,1])
% lgd = legend('$d_1$', '$d_2$', '$d_3$', ...
%      "Location","northeast");
% fontsize(lgd,11,"points");
% grid on;
% set(gcf, "Theme", "light"); % Uncomment for report plots
% 
% Ninit = 5;
% [sig_val, total_events] = bulkETCsim(10,0.01,0.95, Ninit,2, 30, 2);
% hs = T_end./(total_events/Ninit);
% [~,hs_ind] = max(hs);

singleETCsim(0.2, [0;0], 30, 3)

%% Question 4.5

