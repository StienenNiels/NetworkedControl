clear
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
B = [0;1];

A2 = A1/3;
B2 = B1;

despo = [-1-2j, -1+2j];
K1 = place(A1,B1,despo);
K2 = place(A2,B2,despo);

p = zeros(7);
p(1:2,3:5) = [0.49,0.02,0.49;0.49,0.02,0.49];
p(6:7,1:2) = [0.99,0.01;0.99,0.01];
p(3:5,6:7) = [0.99,0.01;0.99,0.01;0.99,0.01];

%% Question 2

h_start = 1e-5;
h_end = 5e-1;
h_steps = 1e2;

h_range = linspace(h_start,h_end,h_steps);

Q = eye(6)*1e-8;

stable = zeros(size(h_range));

i = 1;

% Method 1
warning off
for h = h_range
    [A01, A10, A11, A0j, A1j] = c2d_combi_markov(A1,A2,B,K1,K2,h);
    cvx_begin sdp quiet
        variable P(6,6,7) semidefinite
        subject to
            A0j'*(p(1,3)*P(:,:,3)+p(1,4)*P(:,:,4)+p(1,5)*P(:,:,5))*A0j -P(:,:,1) +Q <= 0;
            A1j'*(p(2,3)*P(:,:,3)+p(2,4)*P(:,:,4)+p(2,5)*P(:,:,5))*A1j -P(:,:,2) +Q <= 0;
            A01'*(p(3,6)*P(:,:,6)+p(3,7)*P(:,:,7)                )*A01 -P(:,:,3) +Q <= 0;
            A11'*(p(4,6)*P(:,:,6)+p(4,7)*P(:,:,7)                )*A11 -P(:,:,4) +Q <= 0;
            A10'*(p(5,6)*P(:,:,6)+p(5,7)*P(:,:,7)                )*A10 -P(:,:,5) +Q <= 0;
            A0j'*(p(6,1)*P(:,:,1)+p(6,2)*P(:,:,2)                )*A0j -P(:,:,6) +Q <= 0;
            A1j'*(p(7,1)*P(:,:,1)+p(7,2)*P(:,:,2)                )*A1j -P(:,:,7) +Q <= 0;
    cvx_end
    stable(i) = strcmp(cvx_status, 'Solved');
    if ~stable(i)
        break
    end
    i = i+1;
end
warning on

% Method 3
warning off
for h = h_range
    [A01, A10, A11, A0j, A1j] = c2d_combi_markov(A1,A2,B,K1,K2,h);
    cvx_begin sdp quiet
        variable P(6,6,7) semidefinite
        subject to
            A01'*p(1,3)*P(:,:,3)*A01 + A11'*p(1,4)*P(:,:,4)*A11 + A10'*p(1,5)*P(:,:,5)*A10 -P(:,:,1) +Q <= 0;
            A01'*p(2,3)*P(:,:,3)*A01 + A11'*p(2,4)*P(:,:,4)*A11 + A10'*p(2,5)*P(:,:,5)*A10 -P(:,:,2) +Q <= 0;

            A0j'*p(3,6)*P(:,:,6)*A0j + A1j'*p(3,7)*P(:,:,7)*A1j -P(:,:,3) +Q <= 0;
            A0j'*p(4,6)*P(:,:,6)*A0j + A1j'*p(4,7)*P(:,:,7)*A1j -P(:,:,4) +Q <= 0;
            A0j'*p(5,6)*P(:,:,6)*A0j + A1j'*p(5,7)*P(:,:,7)*A1j -P(:,:,5) +Q <= 0;

            A0j'*p(6,1)*P(:,:,1)*A0j + A1j'*p(6,2)*P(:,:,2)*A1j -P(:,:,6) +Q <= 0;
            A0j'*p(7,1)*P(:,:,1)*A0j + A1j'*p(7,2)*P(:,:,2)*A1j -P(:,:,7) +Q <= 0;
    cvx_end
    stable(i) = strcmp(cvx_status, 'Solved');
    if ~stable(i)
        break
    end
    i = i+1;
end
warning on

h_max = find(stable==1,1,"last");
%%

figure(212), clf;
hold on;
grid
stairs(h_range, stable, LineWidth=1.5);
plot(h_range(h_max+1), 1,'.', "MarkerSize",20);
lgd = legend("Closed-loop stability", ...
      ['$h_{max} = ', num2str(h_range(h_max), '%.4f'),'$'], "Location","northeast");
fontsize(lgd,12,"points");
ylim([-0.2, 1.2])
xlabel("$h$ [seconds]")
ylabel("Stable (1), Unstable (0)")
set(gcf, "Theme", "light"); % Uncomment for report plots