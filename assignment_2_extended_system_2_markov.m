clear

syms x1 x2 F1 F2

xk = [x1;x2;x2;x2];
q=1;
A11 = F1;
A21 = [0, 1, 0;
       0, 0, 1;
       0, 0, F2];
A22 = [0, 1, 0;
       0, 0, 1;
       0, 0, 1];

A = [A11, [0, 0, 0];
     [0;0;0], A22];
A3 = [A11, [0, 0, 0];
     [0;0;0], A21];
for i = 1:10
if q == 3
    xk = A3*xk
else
    xk = A*xk
end
q = mod(q,3)+1
end