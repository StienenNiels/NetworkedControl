function [A_cl,F,G] = c2d_zoh_FG(A, B, K_static, U_gain, h, tau)
    Fx = expm(A*h);
    G1 = (expm(A*(h-tau)) -eye(2))/A *B;
    Fu = (Fx -eye(2))/A *B -G1;

    F = [Fx, Fu;
        zeros(1,3)];
    G = [G1; 1];
    K = [K_static, U_gain];
    A_cl = F-G*K;
end