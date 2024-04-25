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
K = place(A,B,p);

% % Symbolic definition of ZOH discretization
% syms h
% Fh = expm (A*h);
% Gh = (Fh - eye(2))/ A*B;

% Testing for what ranges of sampling times h the system is stable
h_range = linspace(1e-4,1e0,1e5);
i = 1;
stable = zeros(size(h_range));
lm = zeros(size(h_range));

for h = h_range
    Fh = expm (A*h);
    Gh = (Fh - eye(2))/ A*B;
    lm(i) = max(abs(eig(Fh-Gh*K))); % Spectral radius
    if lm(i) < 1
        stable(i) = 1;
        h_max = h;
    end
    i = i+1;
end
h_index = find(h_range == h_max);

figure(1), clf;
plot(lm,"LineWidth",1.5), hold on;
xline(h_index,"LineWidth",1.5,"Color","k");
legend('Spectral radius',['$h = ', num2str(h_range(h_index), '%.4f'),'$'], "interpreter", "latex");


%% Question 2
% Redo Q1 but now there is a small but constant delay tau in [0,h)
% Lecture 2, slide 12-18

A = [0.3+a-b, 0.5-c;
     0, 1];
B = [0;1];

syms h
Fx = expm (A*h);
Gh = (Fh - eye(2))/ A*B;
