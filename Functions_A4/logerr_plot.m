function logerr_plot(xf_central,xf_compare,filename)
    if nargin < 3
        savefigure = 0;
        plotstate = 1;
    else
        savefigure = 1;
    end

    % Reshape and repeat the centralized solution to be the same size
    norm_xf_cen = norm(xf_central);
    xf_cen = repmat(xf_central, 4, length(xf_compare));
    xf_error = xf_cen-xf_compare;

    % Create a multipaged array to use the pagenorm
    reshaped_matrix = reshape(xf_error, 4, 1, 4, length(xf_compare));
    xf_error_norm = squeeze(pagenorm(reshaped_matrix,2))/norm_xf_cen;
    
    if plotstate
        figure(99),clf
        tiledlayout(1,1)
        nexttile
        hold on
        plot(xf_error_norm')
        yscale('log')
        ylim([1e-3, 1e1])
        % xlim([0,1500])
        legend("$x_{1,f}$", "$x_{2,f}$", "$x_{3,f}$", "$x_{4,f}$")
        xlabel("Iterations")
        ylabel("$\frac{||\textbf{x}_{i,f}-\textbf{x}^*_{f}||_2}{||\textbf{x}^*_{f}||_2}$")
    end

    if savefigure
        figure(99);
        clf;
        tiledlayout(1,1)
        nexttile
        hold on
        plot(xf_error_norm')
        yscale('log')
        ylim([1e-2, 1e1])
        legend("$x_{1,f}$", "$x_{2,f}$", "$x_{3,f}$", "$x_{4,f}$")
        xlabel("Iterations")
        ylabel("$\frac{||\textbf{x}_{i,f}-\textbf{x}^*_{f}||_2}{||\textbf{x}^*_{f}||_2}$")
        grid on
        % set(fig, 'Position', [1000, 100, 500, 400]);
        set(gcf, "Theme", "light");
        filename = "Images\4_" + filename + ".pdf";
        exportgraphics(gcf,filename)
        close(gcf);
    end
end