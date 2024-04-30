clear
close all
clc

%% Initialization of matrices A and B
student_id = 5595738;
a = 5;
b = 9;
c = 8;

A = [0.3+a-b, 0.5-c;
     0, 1];
B = [0;1];

%% Question 1
% Desired poles and pole placement controller
p = [-1-2j, -1+2j];
K_static = place(A,B,p);

% % Symbolic definition of ZOH discretization
% % Lecture 1, slide 41
% syms h
% Fh = expm(A*h);
% Gh = (Fh -eye(2))/A *B;

% Testing for what ranges of sampling times h the system is stable
h_range = linspace(1e-2,1e0,1e5);
i = 1;
stable = zeros(size(h_range));
lm = zeros(size(h_range));

for h = h_range
    Fh = expm(A*h);
    Gh = (Fh -eye(2))/A *B;
    lm(i) = max(abs(eig(Fh-Gh*K_static))); % Spectral radius
    if lm(i) < 1
        stable(i) = 1;
        h_max = h;
    end
    i = i+1;
end

figure(11), clf;
plot(h_range, lm,"LineWidth",1.5), hold on;
plot(h_max,0.999,'.',"MarkerSize",20);
yline(1,"LineWidth",1.5,"Color","w");
xscale("log")
ylim([0.87, 1.3])
lgd = legend('Spectral radius',['$h = ', num2str(h_max, '%.4f'),'$'], "interpreter", "latex", "Location","northwest");
fontsize(lgd,14,"points");
xlabel("$h \; [seconds]$", "Interpreter","latex")
ylabel("$\rho \big(F(h)-G(h)\bar{K}\big)$", "Interpreter","latex")
% set(gcf, "Theme", "light"); % Uncomment for report plots


%% Question 2
% Redo Q1 but now there is a small but constant delay tau in [0,h)
% Lecture 2, slide 12-18

% syms h tau
% Fx = expm(A*h);
% G1 = (expm(A*(h-tau)) -eye(2))/A *B;
% Fu = (Fx -eye(2))/A *B -G1;

% F = [Fx, Fu;
%      zeros(1,3)];
% G = [G1; eye(1)];
% K = [K_static, 0];

% Testing for what ranges of sampling times h the system is stable
h_res = 1e2;
tau_res = 1e2;

h_range = linspace(1e-4,1e0,h_res);
tau_range = linspace(0,5e-1,tau_res);
i = 1;
stable = zeros(size(h_range,2),size(tau_range,2))';
lm = zeros(size(h_range,2),size(tau_range,2))';
h_max = zeros(size(tau_range));

for h = h_range
    j = 1;
    for tau = tau_range
        Fx = expm(A*h);
        G1 = (expm(A*(h-tau)) -eye(2))/A *B;
        Fu = (Fx -eye(2))/A *B -G1;
    
        F = [Fx, Fu;
            zeros(1,3)];
        G = [G1; eye(1)];
        K = [K_static, 0];
        lm(j,i) = max(abs(eig(F-G*K))); % Spectral radius
        if lm(j,i) < 1
            stable(j,i) = 1;
            h_max(i) = h;
        end
        j = j+1;
    end
    i = i+1;
end

% Find the edge of the contour plot
[cline]=contourc(lm,[1 1]);
c_resh = [(1e0-1e-4)/h_res*cline(1,2:end);
          (5e-1)/tau_res*cline(2,2:end)];
ind = find(c_resh(1,:)<c_resh(2,:),1,"first");
poly = polyshape([[0;0], c_resh(:,1:ind)-0.005]');

figure(21), clf;
contour(h_range, tau_range, lm, "ShowText","on", "LineWidth", 1.5, "LabelSpacing", 85), hold on;
plot(poly, "FaceColor", "b", "FaceAlpha", 0.2, "EdgeAlpha",0);
plot(polyshape([0 1 0],[0 1 1]), "FaceColor", "w", "FaceAlpha", 1);
plot([0,1],[0,1], "LineWidth", 1.5, "Color","#D95319");
xlim([0, 1]);
ylim([0, 0.5]);
lgd = legend('Spectral radius contour lines','Stable region','','$\tau = h$', "interpreter", "latex", "Location","northwest");
fontsize(lgd,14,"points");
xlabel("$h \;[seconds]$", "Interpreter","latex")
ylabel("$\tau \;[seconds]$", "Interpreter","latex")
set(gcf, "Theme", "light"); % Uncomment for report plots

%% Question 2.2
% Definitely need to make this plot better
U_gain = 0.9;
K1 = [K_static, 0];
K2 = [K_static, U_gain];

% Testing for what ranges of sampling times h the system is stable
h = 0.4;
tau_range = linspace(0,0.4,1e2);
lm_static = zeros(size(tau_range));
lm_dynamic = zeros(3,size(tau_range,2));

i = 1;
for tau = tau_range
    if tau > h
        break;
    end
    Fx = expm(A*h);
    G1 = (expm(A*(h-tau)) -eye(2))/A *B;
    Fu = (Fx -eye(2))/A *B -G1;

    F = [Fx, Fu;
        zeros(1,3)];
    G = [G1; eye(1)];

    lm_static(i) = max(abs(eig(F-G*K1))); % Spectral radius
    lm_dynamic(:,i) = (abs(eig(F-G*K2))); % Spectral radius
    
    i = i+1;
end

figure(22), clf;
plot(tau_range, lm_static, "LineWidth",1.5), hold on;
plot(tau_range, lm_dynamic, "LineWidth",1.5, "LineStyle","--");
plot(tau_range, max(lm_dynamic), "LineWidth",1.5);
yline(1, "LineWidth",1.5);
% ylim([0.5, 1.4]);
legend('$K=\big[\bar{K} \; 0\big]$','','','',['$K=\big[\bar{K} \;', num2str(U_gain, '%.1f'),'\big]$'],'1', "interpreter", "latex", "Location","northwest");
xlabel("$\tau \;[seconds]$", "Interpreter","latex")
ylabel("$\rho \big(F(h)-G(h)\bar{K}\big)$", "Interpreter","latex")
set(gcf, "Theme", "light"); % Uncomment for report plots


%% Question 3
%% Q3.1
Fx = expm(A*h);
eAds = (expm(A*h) -eye(2))/A *B;
Fu1 = -A\expm(A*h)*B + eAds + A\eAds/h;
Fu2 = A\expm(A*h)*B - A\eAds/h;

F1 = [Fx, Fu1, Fu2; zeros(1,4);0,0,1,0];
G1 = [0;0;1;0];

K = [K_static,0,0];

%% Q3.2
% Testing for what ranges of sampling times h the system is stable
h_res = 1e4;
tau_res = 2e0;

h_range = linspace(1e-4,1e0,h_res);
tau_range = linspace(0,5e-1,tau_res);
i = 1;
stable = zeros(size(h_range,2),size(tau_range,2))';
lm = zeros(size(h_range,2),size(tau_range,2))';
h_max = zeros(size(tau_range));

for h = h_range
    j = 1;
    for tau = tau_range
        Fx = expm(A*h);
        eAds = (expm(A*h) -eye(2))/A *B;
        Fu1 = -A\expm(A*h)*B + eAds + A\eAds/h;
        Fu2 = A\expm(A*h)*B - A\eAds/h;
        
        F1 = [Fx, Fu1, Fu2; zeros(1,4);0,0,1,0];
        G1 = [0;0;1;0];
        
        K = [K_static,0,0];
        lm(j,i) = max(abs(eig(F1-G1*K))); % Spectral radius
        if lm(j,i) < 1
            stable(j,i) = 1;
            h_max(i) = h;
        end
        j = j+1;
    end
    i = i+1;
end

x32 = 0.12426;
poly32 = polyshape([[x32, x32; x32, 0], c_resh(:,1:ind)-0.005]');

figure(32), clf;
contour(h_range, tau_range, lm, "ShowText","on", "LineWidth", 1.5, "LabelSpacing", 85), hold on;
plot(poly32, "FaceColor", "b", "FaceAlpha", 0.2, "EdgeAlpha",0);
plot(polyshape([0 x32 x32],[0 x32 0]), "FaceColor", "g", "FaceAlpha", 0.45, "EdgeAlpha",0);
plot(polyshape([0 1 0],[0 1 1]), "FaceColor", "w", "FaceAlpha", 1);
plot([0,1],[0,1], "LineWidth", 1.5, "Color","#D95319");
xlim([0, 1]);
ylim([0, 0.5]);
lgd = legend('Spectral radius contour lines','Stable region $Q_2$','Stable region $Q_3$','','$\tau = h$', "interpreter", "latex", "Location","northwest");
fontsize(lgd,14,"points");
xlabel("$h \;[seconds]$", "Interpreter","latex")
ylabel("$\tau \;[seconds]$", "Interpreter","latex")
set(gcf, "Theme", "light"); % Uncomment for report plots

%% Q3.3
U1_gain = 0.8;
U2_gain = 0.5;

h_res = 1e4;
tau_res = 2;

h_range = linspace(1e-4,1e0,h_res);
tau_range = linspace(0,1,tau_res);
i = 1;
stable = zeros(size(h_range,2),size(tau_range,2))';
lm = zeros(size(h_range,2),size(tau_range,2))';
h_max = zeros(size(tau_range));

for h = h_range
    j = 1;
    for tau = tau_range
        Fx = expm(A*h);
        eAds = (expm(A*h) -eye(2))/A *B;
        Fu1 = -A\expm(A*h)*B + eAds + A\eAds/h;
        Fu2 = A\expm(A*h)*B - A\eAds/h;
        
        F1 = [Fx, Fu1, Fu2; zeros(1,4);0,0,1,0];
        G1 = [0;0;1;0];
        
        K = [K_static,U1_gain,U2_gain];
        lm(j,i) = max(abs(eig(F1-G1*K))); % Spectral radius
        if lm(j,i) < 1
            stable(j,i) = 1;
            h_max(i) = h;
        end
        j = j+1;
    end
    i = i+1;
end

idx = find(lm(1,:)<1,1,"last");
x33 = h_range(idx);

figure(33), clf;
contour(h_range, tau_range, lm, "ShowText","on", "LineWidth", 1.5, "LabelSpacing", 85, "LevelList",[1 1.3 1.8]), hold on;
plot(polyshape([0 x32 x32],[0 x32 0]), "FaceColor", "g", "FaceAlpha", 0.45, "EdgeAlpha",0);
plot(polyshape([x32 x32 x33 x33],[0 x32 x33 0]), "FaceColor", "b", "FaceAlpha", 0.2, "EdgeAlpha",0);
plot(polyshape([0 1 0],[0 1 1]), "FaceColor", "w", "FaceAlpha", 1);
plot([0,1],[0,1], "LineWidth", 1.5, "Color","#D95319");
xlim([0, 1]);
ylim([0, 1]);
lgd = legend('Spectral radius contour lines','Stable region $K=\big[\bar{K} \; 0.0\; 0.0\big]$',['Stable region $K=\big[\bar{K} \;', num2str(U1_gain, '%.1f'),'\;',num2str(U2_gain, '%.1f'),'\big]$'],'','$\tau = h$', "interpreter", "latex", "Location","northwest");
fontsize(lgd,14,"points");
xlabel("$h \;[seconds]$", "Interpreter","latex")
ylabel("$\tau \;[seconds]$", "Interpreter","latex")
set(gcf, "Theme", "light"); % Uncomment for report plots

%% Question 4
% Find a common Lyapunov function that holds for both zoh and linear order
% interpolation

%% Q4.2 % Check met Alex

% Testing for what ranges of sampling times h the system is stable
h_res = 1e3;
tau_res = 1e3;

h_range = linspace(1e-4,1e0,h_res);
tau_range = linspace(0,1,tau_res);
i = 1;
stable = zeros(size(h_range,2),size(tau_range,2))';
lmzoh = zeros(size(h_range,2),size(tau_range,2))';
lmfoh = zeros(size(h_range,2),size(tau_range,2))';
lmsq1 = zeros(size(h_range,2),size(tau_range,2))';
lmsq2 = zeros(size(h_range,2),size(tau_range,2))';
h_max = zeros(size(tau_range));

for h = h_range
    j = 1;
    for tau = tau_range
        if tau > h
            break;
        end
        Fx = expm(A*h);
        G1 = (expm(A*(h-tau)) -eye(2))/A *B;
        Fu = (Fx -eye(2))/A *B -G1;
    
        F = [Fx, Fu, [0;0];
            zeros(2,4)];
        G = [G1; 1;0];
        K = [K_static, U_gain,0];
        lmzoh(j,i) = max(abs(eig(F-G*K))); % Spectral radius

        Fx = expm(A*h);
        eAds = (expm(A*h) -eye(2))/A *B;
        Fu1 = -A\expm(A*h)*B + eAds + A\eAds/h;
        Fu2 = A\expm(A*h)*B - A\eAds/h;
        
        F2 = [Fx, Fu1, Fu2; zeros(1,4);0,0,1,0];
        G2 = [0;0;1;0];
        
        K2 = [K_static,U1_gain,U2_gain];
        lmfoh(j,i) = max(abs(eig(F2-G2*K2))); % Spectral radius
        lmsq1(j,i) = max(abs(eig( (F-G*K)*(F2-G2*K2) )));
        lmsq2(j,i) = max(abs(eig( (F-G*K)*(F2-G2*K2)*(F2-G2*K2) )));
        j = j+1;
    end
    i = i+1;
end

figure(431), clf;
% contour(h_range, tau_range, lmzoh, "ShowText","on", "LineWidth", 1.5, "LabelSpacing", 85), hold on;
% contour(h_range, tau_range, lmfoh, "ShowText","on", "LineWidth", 1.5, "LabelSpacing", 85), hold on;
contour(h_range, tau_range, lmsq1, "ShowText","on", "LineWidth", 1.5, "LabelSpacing", 85), hold on;
% contour(h_range, tau_range, lmsq2, "ShowText","on", "LineWidth", 1.5, "LabelSpacing", 85), hold on;
% plot(polyshape([0 x32 x32],[0 x32 0]), "FaceColor", "g", "FaceAlpha", 0.45, "EdgeAlpha",0);
% plot(polyshape([x32 x32 x33 x33],[0 x32 x33 0]), "FaceColor", "b", "FaceAlpha", 0.2, "EdgeAlpha",0);
% plot(polyshape([0 1 0],[0 1 1]), "FaceColor", "w", "FaceAlpha", 1);
plot([0,1],[0,1], "LineWidth", 1.5, "Color","#D95319");
xlim([0, 1]);
ylim([0, 1]);
lgd = legend('sq1','$\tau = h$', "interpreter", "latex", "Location","northwest");
fontsize(lgd,14,"points");
xlabel("$h \;[seconds]$", "Interpreter","latex")
ylabel("$\tau \;[seconds]$", "Interpreter","latex")
% set(gcf, "Theme", "light"); % Uncomment for report plots

figure(432), clf;
% contour(h_range, tau_range, lmzoh, "ShowText","on", "LineWidth", 1.5, "LabelSpacing", 85), hold on;
% contour(h_range, tau_range, lmfoh, "ShowText","on", "LineWidth", 1.5, "LabelSpacing", 85), hold on;
% contour(h_range, tau_range, lmsq1, "ShowText","on", "LineWidth", 1.5, "LabelSpacing", 85), hold on;
contour(h_range, tau_range, lmsq2, "ShowText","on", "LineWidth", 1.5, "LabelSpacing", 85), hold on;
% plot(polyshape([0 x32 x32],[0 x32 0]), "FaceColor", "g", "FaceAlpha", 0.45, "EdgeAlpha",0);
% plot(polyshape([x32 x32 x33 x33],[0 x32 x33 0]), "FaceColor", "b", "FaceAlpha", 0.2, "EdgeAlpha",0);
% plot(polyshape([0 1 0],[0 1 1]), "FaceColor", "w", "FaceAlpha", 1);
plot([0,1],[0,1], "LineWidth", 1.5, "Color","#D95319");
xlim([0, 1]);
ylim([0, 1]);
lgd = legend('sq2','$\tau = h$', "interpreter", "latex", "Location","northwest");
fontsize(lgd,14,"points");
xlabel("$h \;[seconds]$", "Interpreter","latex")
ylabel("$\tau \;[seconds]$", "Interpreter","latex")
% set(gcf, "Theme", "light"); % Uncomment for report plots

%% Question 5

% Constant h and no delay
% Initilization of systems:
h1 = 0.4;
A1 = [0.3+a-b, 0.5-c;
     0, 1];
B1 = [0;1];

h2 = 3*h1;
A2 = A1/3;
B2 = B1;

% Packets drop out alternating
% So system 1 drops a packet every 6th sample
% And system 2 drops every 2nd packet
% First attempt simply use the same pole placement for both
p = [-1-2j, -1+2j];
K_1 = place(A1,B1,p);
K_2 = place(A2,B2,p);

% to-zero analysis



% to-hold analysis




function A_cl = c2d_zoh(A,B,K_static,U_gain,h,tau)
    Fx = expm(A*h);
    G1 = (expm(A*(h-tau)) -eye(2))/A *B;
    Fu = (Fx -eye(2))/A *B -G1;

    F = [Fx, Fu;
        zeros(1,3)];
    G = [G1; eye(1)];
    K = [K_static, U_gain];
    A_cl = F-G*K;
end