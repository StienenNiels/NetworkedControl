clear
clc

addpath("Functions_A4\")
addpath("Images\")

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% state vector is (x,y,xdot,ydot)
aircraft; % Load the given parameters
Planes = planes_gen(); % Define structures

% Create structures for each plane to make workspace more structured
for i = 1:4
    Planes(i).plane = i;
    Planes(i).Tf = Tfinal;
    Planes(i).umax = umax;
    Planes(i).A = eval(sprintf('A%d;', i));
    Planes(i).B = eval(sprintf('B%d;', i));
    Planes(i).x0 = eval(sprintf('x0%d;', i));
    Planes(i) = predmodgen(Planes(i));
    Planes(i) = optgen(Planes(i));
end
dim = Planes(1).dim;
clearvars -except Planes Tfinal umax dim

%% Centralized solution
central_sol(Planes, 1);