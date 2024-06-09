function xf_error_norm = err_norm(xf_central,xf_compare,average)

    % Reshape and repeat the centralized solution to be the same size
    norm_xf_cen = norm(xf_central);
    xf_cen = repmat(xf_central, 4, length(xf_compare));
    xf_error = xf_cen-xf_compare;

    % Create a multipaged array to use the pagenorm
    reshaped_matrix = reshape(xf_error, 4, 1, 4, length(xf_compare));
    xf_error_norm = squeeze(pagenorm(reshaped_matrix,2))/norm_xf_cen;
    
    if average
        xf_error_norm = mean(xf_error_norm);
    end
end