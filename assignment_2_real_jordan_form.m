clear


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
U_gain = 0.9;

[Q,J] = jordan(A)
S1 = [1,0;0,0];
S2 = [0,0;0,1];

lambda_1 = J(1,1);
lambda_2 = J(2,2);

W1 = Q\S1*Q;
W2 = Q\S2*Q;
% eAs = exp(J(1,1)*s)*W1 + exp(J(2,2)*s)*W2
eAs = exp(lambda_1*s)*W1 + exp(lambda_2*s)*W2
F = [subs(eAs,s,h), int(eAs,s,ht,h)*B; 0,0,0]
G = [int(eAs,s,0,ht)*B;1];

alpha1 = exp(ht*lambda_1)
alpha2 = exp(ht*lambda_2)

F0 = subs(F,[alpha1,alpha2],[0,0])
simplify(F-F0);
F1 = diff(subs(simplify(F),alpha1,C),C)
simplify(F-F0-F1*alpha1);
F2 = diff(subs(simplify(F),alpha2,C),C)
simplify(F-F0-F1*alpha1-F2*alpha2);

% Check 
% F0 + alpha_1*F1 + alpha_2*F2
simplify(F - (F0 + alpha1*F1 + alpha2*F2))
tau_range = 0:0.01:1;
i=1;

vals = zeros(2,size(tau_range,2));

for tau = tau_range
    vals(:,i) = [exp(-3.7*(1-tau));exp(1*(1-tau))];
    i = i+1;
end


figure(1),clf;
plot(vals(1,:),vals(2,:));
xlabel("$\alpha_1(\tau)$")
ylabel("$\alpha_2(\tau)$")

% % We have 2 eigenvalues in our system so 2 alpha's
% 
% alpha_1(tau_k) = ((h-tau_k))*exp(J(1,1)*(h-tau_k))
% alpha_2(tau_k) = ((h-tau_k))*exp(J(2,2)*(h-tau_k))
