function rho_opt_plot(xf,varstruc,legendtype,filename)
    if nargin < 4
        savefigure = 0;
        plotstate = 1;
    else
        plotstate = 1;
        savefigure = 1;
    end
    
    fields = fieldnames(varstruc);
    vars = varstruc.(fields{1});
    nvars = length(varstruc.(fields{1}));

    rho_length = zeros(1,nvars);
    for i = 1:nvars
        rho_length(i) = length(xf{i});
    end
    rho_length = filloutliers(rho_length,"linear","movmedian",5);
    rho_length = sgolayfilt(rho_length,3,33);


    [~, ind] = min(rho_length);

    if plotstate
        % existing_figures = findobj('Type', 'figure');
        % if isempty(existing_figures)
        %     figure(1);
        % else
        %     figure_numbers = arrayfun(@(x) x.Number, existing_figures);
        %     new_fig_number = max(figure_numbers) + 1;
        %     figure(new_fig_number);
        % end
        figure(99);
        clf
        tiledlayout(1,1)
        nexttile
        hold on
        plot(vars, rho_length')
        xline(6.5, "-.", "LineWidth",1.5)
        xline(vars(ind), "-.", "LineWidth",1.5)
        ylim([0,200])
        xlim([0,max(vars)])
        legend("", "$\rho = 6.5$", "$\rho = "+vars(ind)+"$")
        xlabel("$\rho$")
        ylabel("Iterations")
        title("Iterations needed for $\Delta x_f \leq 10^{-5}$ for varying $\rho$")
        grid on
        % set(gcf, 'Position', get(0, 'Screensize'));
    end

    if savefigure
        set(gcf, 'Position', [1000, 100, 500, 400]);
        set(gcf, "Theme", "light");
        filename = "Images\4_" + filename + ".pdf";
        exportgraphics(gcf,filename)
        close(gcf);
    end
end