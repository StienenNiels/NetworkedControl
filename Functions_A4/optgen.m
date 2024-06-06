function [plane]=optgen(plane)

dim.nx = size(plane.A,1);
dim.nu = size(plane.B,2);
dim.N  = plane.Tf;
S = plane.S;
T = plane.T;
x0 = plane.x0;
umax = plane.umax;
Tf = plane.Tf;

% Cost function
% u_opt = u^T*H*u + 2*h'*u
plane.H = S'*S + eye(dim.N*dim.nu);
plane.h = (S'*T*x0);

% Plane final state equality
plane.A_eq = S(dim.nx*(dim.N-1)+1:end,:);
plane.b_eq = T(dim.nx*(dim.N-1)+1:end,:)*x0;

% Control input constraint inequality
plane.A_u = kron(eye(dim.N*dim.nu), [-1;1]);
plane.b_u = repmat([umax/Tf; umax/Tf],dim.N*dim.nu,1);
 
end
