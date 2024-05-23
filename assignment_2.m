clear
close all
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

A1 = [0.3+a-b, 0.5-c;
     0, 1];
B1 = [0;1];

A2 = A1/3;
B2 = B1;

p = [-1-2j, -1+2j];
K1 = place(A1,B1,p);
K2 = place(A2,B2,p);

%% Question 1

h_start = 1e-5;
h_end = 5e-1;
h_steps = 1e2;

h_range = linspace(h_start,h_end,h_steps);

Qh = eye(3);
Qz = eye(2);

stableh = zeros(size(h_range));
stablez = zeros(size(h_range));

i = 1;

for h = h_range
    [A_zp, A_znp, A_hp, A_hnp] = c2d_zero_hold(A1, B1, K1, h);
    Ah1 = A_hp*A_hp*A_hp;
    Ah2 = A_hp*A_hp*A_hnp;
    cvx_begin sdp quiet
        variable Ph(3,3) semidefinite
        subject to
            Ah1'*Ph*Ah1 -Ph +Qh <= 0;
            Ah2'*Ph*Ah2 -Ph +Qh <= 0;
    cvx_end
    stableh(i) = strcmp(cvx_status, 'Solved'); %any(isnan(Ph), 'all');
    Az1 = A_zp*A_zp*A_zp;
    Az2 = A_zp*A_zp*A_znp;
    cvx_begin sdp quiet
        variable Pz(2,2) semidefinite
        subject to
            Az1'*Pz*Az1 -Pz +Qz <= 0;
            Az2'*Pz*Az2 -Pz +Qz <= 0;
    cvx_end
    stablez(i) = strcmp(cvx_status, 'Solved'); %any(isnan(Pz), 'all');
    i = i+1;
end

h_max_hold = find(stableh==1,1,"last");
h_max_zero = find(stablez==1,1,"last");

%%
figure(1), clf;
hold on;
grid
stairs(h_range, stableh, LineWidth=1.5);
stairs(h_range, stablez, "--", LineWidth=1.5);
plot(h_range(h_max_hold+1), 1,'.', "MarkerSize",20);
plot(h_range(h_max_zero+1), 1,'.', "MarkerSize",20);
lgd = legend("Closed-loop stability to-hold", "Closed-loop stability to-zero", ...
      ['$h_{max} = ', num2str(h_range(h_max_hold), '%.4f'),'$'], ['$h_{max} = ', num2str(h_range(h_max_zero), '%.4f'),'$'],"Location","southwest");
fontsize(lgd,11,"points");
ylim([-0.2, 1.2])
xlabel("$h$ [seconds]")
ylabel("Stable (1), Unstable (0)")
set(gcf, "Theme", "light"); % Uncomment for report plots