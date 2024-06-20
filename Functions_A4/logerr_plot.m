function logerr_plot(xf_central,xf_compare,filename)
    if nargin < 3
        savefigure = 0;
        plotstate = 1;
    else
        plotstate = 1;
        savefigure = 1;
    end

    xf_error_norm = err_norm(xf_central,xf_compare,0);
    
    if plotstate
        figure(99),clf
        tiledlayout(1,1)
        nexttile
        hold on
        plot(xf_error_norm')
        yscale('log')
        ylim([1e-10, 1e1])
        % xlim([0,1500])
        legend("$x_{1,f}$", "$x_{2,f}$", "$x_{3,f}$", "$x_{4,f}$")
        xlabel("Iterations")
        ylabel("$\frac{||\textbf{x}_{i,f}-\textbf{x}^*_{f}||_2}{||\textbf{x}^*_{f}||_2}$")
        title("Logarithmic error sequence")
        grid on
    end

    if savefigure
        set(gcf, 'Position', [1000, 100, 500, 400]);
        set(gcf, "Theme", "light");
        filename = "Images\4_" + filename + ".pdf";
        exportgraphics(gcf,filename)
        close(gcf);
    end
end