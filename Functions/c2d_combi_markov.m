function [A01, A10, A11, A0j, A1j] = c2d_combi_markov(A1,A2,B,K1,K2,h)

    [~, ~, A1_hp, A1_hnp] = c2d_zero_hold(A1, B, K1, h);
    [~, ~, A2_hp, A2_hnp] = c2d_zero_hold(A2, B, K2, 3*h);

    zero = zeros(size(A1_hp));
    A01 = [A1_hp, zero; zero, A2_hnp];
    A10 = [A1_hnp, zero; zero, A2_hp];
    A11 = [A1_hnp, zero; zero, A2_hnp];
    A0j = [A1_hp, zero; zero, eye(size(A2_hp))];
    A1j = [A1_hnp, zero; zero, eye(size(A2_hp))];

end