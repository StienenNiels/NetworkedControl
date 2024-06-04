clear
clc

addpath("Functions\")
addpath("Images\")

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% state vector is (x,y,xdot,ydot)
aircraft; % Load the given parameters