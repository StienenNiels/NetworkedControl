clear
clc

addpath("Functions\")

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

syms h tau_k i lambda_i lambda_1 lambda_2 s ht C alpha_1 alpha_2
student_id = 5595738;
a = 5;
b = 9;
c = 8;

A = [0.3+a-b, 0.5-c;
     0, 1];
B = [0;1];

p = [-1-2j, -1+2j];
K_static = place(A,B,p);
K = [K_static,0];
U_gain = 0.9;

[Q,J] = jordan(A);
sym(Q);
sym(J);
S1 = [1,0;0,0];
S2 = [0,0;0,1];

lambda_1 = J(1,1);
lambda_2 = J(2,2);

W1 = Q*S1/Q;
W2 = Q*S2/Q;
% eAs = exp(J(1,1)*s)*W1 + exp(J(2,2)*s)*W2
eAs = exp(lambda_1*s)*W1 + exp(lambda_2*s)*W2;
F = [subs(eAs,s,h), int(eAs,s,ht,h)*B; 0,0,0];
G = [int(eAs,s,0,ht)*B;1];

alpha1 = exp(ht*lambda_1);
alpha2 = exp(ht*lambda_2);

F0 = subs(F,[alpha1,alpha2],[0,0]);
simplify(F-F0);
F1 = diff(subs(simplify(F),alpha1,C),C);
simplify(F-F0-F1*alpha1);
F2 = diff(subs(simplify(F),alpha2,C),C);
simplify(F-F0-F1*alpha1-F2*alpha2);

G0 = subs(G,[alpha1,alpha2],[0,0]);
simplify(G-G0);
G1 = diff(subs(simplify(G),alpha1,C),C);
simplify(G-G0-G1*alpha1);
G2 = diff(subs(simplify(G),alpha2,C),C);
simplify(G-G0-G1*alpha1-G2*alpha2);

% Check 
% F0 + alpha_1*F1 + alpha_2*F2
simplify(F - (F0 + alpha1*F1 + alpha2*F2));
tau_range = 0:0.01:1;
i=1;

vals = zeros(2,size(tau_range,2));

for tau = tau_range
    vals(:,i) = [exp(-3.7*(1-tau));exp(1*(1-tau))];
    i = i+1;
end

poly = polyshape([min(vals(1,:)) min(vals(1,:)) max(vals(1,:))],[max(vals(2,:)) min(vals(2,:)) min(vals(2,:))]);
figure(1),clf;
hold on
grid on;
plot(vals(1,:),vals(2,:), LineWidth=1.5);
plot(poly, "FaceColor", "b", "FaceAlpha", 0.2, "EdgeAlpha",0);
xlabel("$\alpha_1(\tau)$")
ylabel("$\alpha_2(\tau)$")
lgd = legend("$\alpha_1$ vs $\alpha_2$", "Polytopic overapproximation");
fontsize(lgd,14,"points");
set(gcf, "Theme", "light"); % Uncomment for report plots

%% Question 3.2
h_start = 0.15;
h_end = 0.2;
h_steps = 1e2;

h_range = linspace(h_start,h_end,h_steps);

Qh = eye(3);
Qz = eye(2);

stable = zeros(size(h_range));

i = 1;
syms h tau alpha1 alpha2
Fprox(h,alpha1, alpha2) = F0 + alpha1*F1 + alpha2*F2;
Gprox(h,alpha1, alpha2) = G0 + alpha2*G1 + alpha2*G2;

LMI_set = matlabFunction(Fprox, Gprox, File="LMI_finder", Vars=[h,alpha1,alpha2]);

for h = h_range
    a1max = exp(-3.7*eps(0));
    a1min = exp(-3.7*h);
    a2min = exp(eps(0));
    a2max = exp(h);
    [Fp1,Gp1] = LMI_set(h,a1max,a2min);
    [Fp2,Gp2] = LMI_set(h,a1min,a2min);
    [Fp3,Gp3] = LMI_set(h,a1min,a2max);
    Ah1 = Fp1-Gp1*K;
    Ah2 = Fp2-Gp2*K;
    Ah3 = Fp3-Gp3*K;
    cvx_begin sdp quiet
        variable Ph(3,3) semidefinite
        subject to
            Ah1'*Ph*Ah1 -Ph +Qh <= 0;
            Ah2'*Ph*Ah2 -Ph +Qh <= 0;
            Ah3'*Ph*Ah3 -Ph +Qh <= 0;
    cvx_end
    stable(i) = strcmp(cvx_status, 'Solved'); %any(isnan(Ph), 'all');
    if stable(i) == 0
        break;
    end
    i = i+1;
end

h_max_ind = find(stable==1,1,"last");
h_max = h_range(h_max_ind)

poly = polyshape([0, h_max, h_max],[0,0,h_max]);

figure(21)
hold on;
plot(poly, "FaceColor", "b", "FaceAlpha", 0.5, "EdgeAlpha",0);
xlim([0, 1]);
ylim([0, 0.5]);
lgd = legend('Spectral radius contour lines','Stable region','','$\tau = h$','Polytopic over-approximation result', "Location","northwest");
fontsize(lgd,14,"points");
xlabel("$h \;[seconds]$", "Interpreter","latex")
ylabel("$\tau \;[seconds]$", "Interpreter","latex")
set(gcf, "Theme", "light"); % Uncomment for report plots

%% Question 3.3

slope = -(max(vals(2,:))-min(vals(2,:))) /(max(vals(1,:))-min(vals(1,:)));
 
offset = 1.89;
syms x
y0 = double(solve(slope*x+offset==1,x));

poly = polyshape([min(vals(1,:)) min(vals(1,:)) y0 max(vals(1,:))],[max(vals(2,:)) (offset+slope*min(vals(1,:))) min(vals(2,:)) min(vals(2,:))]);
figure(33),clf;
hold on
grid on;
plot(vals(1,:),vals(2,:), LineWidth=1.5);
fplot(slope*x+offset,[min(vals(1,:)) y0], LineWidth=1.5);
plot(poly, "FaceColor", "b", "FaceAlpha", 0.2, "EdgeAlpha",0);
xlabel("$\alpha_1(\tau)$")
ylabel("$\alpha_2(\tau)$")
lgd = legend("$\alpha_1$ vs $\alpha_2$","Tangent line", "Polytopic overapproximation");
fontsize(lgd,14,"points");
set(gcf, "Theme", "light"); % Uncomment for report plots

% stable = zeros(size(h_range));
% i = 1;
% 
% for h = h_range
%     a1max = exp(-3.7*eps(0));
%     a1cep = y0;
%     a1min = exp(-3.7*h);
%     a2min = exp(eps(0));
%     a2cep = (offset+slope*min(vals(1,:)));
%     a2max = exp(h);
% 
%     [Fp1,Gp1] = LMI_set(h,a1max,a2min);
%     [Fp2,Gp2] = LMI_set(h,a1min,a2min);
%     [Fp3,Gp3] = LMI_set(h,a1min,a2max);
%     [Fp4,Gp4] = LMI_set(h,a1min,a2max);
%     Ah1 = Fp1-Gp1*K;
%     Ah2 = Fp2-Gp2*K;
%     Ah3 = Fp3-Gp3*K;
%     Ah4 = Fp4-Gp4*K;
%     cvx_begin sdp quiet
%         variable Ph(3,3) semidefinite
%         subject to
%             Ah1'*Ph*Ah1 -Ph +Qh <= 0;
%             Ah2'*Ph*Ah2 -Ph +Qh <= 0;
%             Ah3'*Ph*Ah3 -Ph +Qh <= 0;
%             Ah4'*Ph*Ah4 -Ph +Qh <= 0;
%     cvx_end
%     stable(i) = strcmp(cvx_status, 'Solved'); %any(isnan(Ph), 'all');
%     if stable(i) == 0
%         break;
%     end
%     i = i+1;
% end
% 
% h_max_ind = find(stable==1,1,"last");
% h_max = h_range(h_max_ind)
% 
% poly = polyshape([0, h_max, h_max],[0,0,h_max]);
% 
% figure(21)
% hold on;
% plot(poly, "FaceColor", "b", "FaceAlpha", 0.5, "EdgeAlpha",0);
% xlim([0, 1]);
% ylim([0, 0.5]);
% lgd = legend('Spectral radius contour lines','Stable region','','$\tau = h$','Polytopic over-approximation result', "Location","northwest");
% fontsize(lgd,14,"points");
% xlabel("$h \;[seconds]$", "Interpreter","latex")
% ylabel("$\tau \;[seconds]$", "Interpreter","latex")
% set(gcf, "Theme", "light"); % Uncomment for report plots
