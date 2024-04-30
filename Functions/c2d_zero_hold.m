function [A_zp, A_znp, A_hp, A_hnp] = c2d_zero_hold(A, B, K, h)
  
F = expm(A*h);
G = (expm(A*h) -eye(2))/A *B;

% to zero closed loop dynamics
A_zp = F-G*K;
A_znp = F;

% to hold closed loop dynamics
A_hp = [F-G*K, [0;0];
        -K, 0];
A_hnp = [F, G;
         0,0,1];
end