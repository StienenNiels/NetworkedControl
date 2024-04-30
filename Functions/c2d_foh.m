function A_cl = c2d_foh(A,B,K_static,U1_gain, U2_gain, h, tau)
    Fx = expm(A*h);
    eAds = (expm(A*h) -eye(2))/A *B;
    Fu1 = -A\B + A\eAds/h;
    Fu2 = A\expm(A*h)*B - A\eAds/h;
    
    F1 = [Fx, Fu1, Fu2; zeros(1,4);0,0,1,0];
    G1 = [0;0;1;0];
    
    K = [K_static,U1_gain,U2_gain];

    A_cl = F1-G1*K;
end