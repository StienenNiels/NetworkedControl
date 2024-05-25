function [S00, S01, S10, S11] = c2d_combi(A1,A2,B,K1,K2,h,tozero)

    [A1_zp, A1_znp, A1_hp, A1_hnp] = c2d_zero_hold(A1, B, K1, h);
    [A2_zp, A2_znp, A2_hp, A2_hnp] = c2d_zero_hold(A2, B, K2, 3*h);
    
    if tozero
        zero = zeros(size(A1_zp));
        A_00 = [A1_zp, zero; zero, A2_zp];
        A_01 = [A1_zp, zero; zero, A2_znp];
        A_10 = [A1_znp, zero; zero, A2_zp];
        A_11 = [A1_znp, zero; zero, A2_znp];
        A_nh = [A1_zp, zero; zero, eye(size(A2_zp))];
    else
        zero = zeros(size(A1_hp));
        A_00 = [A1_hp, zero; zero, A2_hp];
        A_01 = [A1_hp, zero; zero, A2_hnp];
        A_10 = [A1_hnp, zero; zero, A2_hp];
        A_11 = [A1_hnp, zero; zero, A2_hnp];
        A_nh = [A1_hp, zero; zero, eye(size(A2_hp))];
    end
    
    S00 = A_00*A_nh*A_nh;
    S01 = A_01*A_nh*A_nh;
    S10 = A_10*A_nh*A_nh;
    S11 = A_11*A_nh*A_nh;

end