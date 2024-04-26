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
K = [K_static, 0.5];

% Testing for what ranges of sampling times h the system is stable
h_res = 1e3;
tau_res = 1e3;

h_range = 0.4;
tau_range = linspace(0,5e-1,tau_res);
i = 1;
stable = zeros(size(h_range,2),size(tau_range,2))';
lm = zeros(size(h_range,2),size(tau_range,2))';
h_max = zeros(size(tau_range));

j = 1;
for tau = tau_range
    if tau >= h
        break;
    end
    Fx = expm(A*h);
    G1 = (expm(A*(h-tau)) -eye(2))/A *B;
    Fu = (Fx -eye(2))/A *B -G1;

    F = [Fx, Fu;
        zeros(1,3)];
    G = [G1; eye(1)];
    K = [K_static, 0.5];
    lm(j,i) = max(abs(eig(F-G*K))); % Spectral radius
    if lm(j,i) < 1
        stable(j,i) = 1;
        h_max(i) = h;
    end
    j = j+1;
end
% h_index = find(h_range == h_max);

figure(22), clf;
plot(tau_range, lm);
legend('Stable region boundary', "interpreter", "latex");
xlabel("$\tau [seconds]$", "Interpreter","latex")
ylabel("$\rho$", "Interpreter","latex")
yline(1);

%% Question 3

% Fx = expm(A*h);
% eAds = (expm(A*h) -eye(2))/A *B;
% Fu1 = -A\expm(A*h)*B + eAds + A\eAds/h;
% Fu2 = A\expm(A*h)*B - A\eAds/h;
% 
% F1 = [Fx, Fu1, Fu2; zeros(1,4);0,0,1,0];
% G1 = [0;0;1;0];
% 
% K = [K_static,0,0];

% Testing for what ranges of sampling times h the system is stable
h_range = linspace(1e-4,1e0,1e5);
i = 1;
stable = zeros(size(h_range));
lm = zeros(size(h_range));

for h = h_range
    Fx = expm(A*h);
    eAds = (expm(A*h) -eye(2))/A *B;
    Fu1 = -A\expm(A*h)*B + eAds + A\eAds/h;
    Fu2 = A\expm(A*h)*B - A\eAds/h;
    
    F1 = [Fx, Fu1, Fu2; zeros(1,4);0,0,1,0];
    G1 = [0;0;1;0];
    
    K = [K_static,0,0];
    lm(i) = max(abs(eig(F1-G1*K))); % Spectral radius
    if lm(i) < 1
        stable(i) = 1;
        h_max = h;
    end
    i = i+1;
end
h_index = find(h_range == h_max);

figure(32), clf;
plot(lm,"LineWidth",1.5), hold on;
xline(h_index,"LineWidth",1.5,"Color","k");
legend('Spectral radius',['$h = ', num2str(h_range(h_index), '%.4f'),'$'], "interpreter", "latex");


%% Question 4
% Find a common Lyapunov function that holds for both zoh and linear order
% interpolation