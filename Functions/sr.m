function spectral_radius = sr(A_cl)
% Calculates the spectral radius
spectral_radius = max(abs(eig(A_cl)));
end