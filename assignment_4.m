clear
clc

addpath("Functions\")
addpath("Images\")

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% state vector is (x,y,xdot,ydot)
aircraft; % Load the given parameters

% Create structures for each plane to make workspace more structured
for i = 1:445
    eval(sprintf('P%d.A = A%d;', i, i));
    eval(sprintf('P%d.B = B%d;', i, i));
    eval(sprintf('P%d.x0 = x0%d;', i, i));
end
clearvars -except P1 P2 P3 P4 Tfinal umax


% Objective
% Write function to generate prediction matrices for Np steps