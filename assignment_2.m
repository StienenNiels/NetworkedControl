clear
close all
clc

addpath("Functions\")

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
