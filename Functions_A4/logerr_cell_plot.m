function logerr_cell_plot(xf,varstruc,legendtype,filename)
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
    if strcmp(legendtype, "nesterov") || strcmp(legendtype, "dual_comp") ...
            || strcmp(legendtype, "consensus")
        vars2 = varstruc.(fields{2});
    end

    if plotstate
        existing_figures = findobj('Type', 'figure');
        if isempty(existing_figures)
            figure(1);
        else
            figure_numbers = arrayfun(@(x) x.Number, existing_figures);
            new_fig_number = max(figure_numbers) + 1;
            figure(new_fig_number);
        end
        clf
        tiledlayout(1,1)
        nexttile
        hold on
        for i = 1:nvars
            if strcmp(legendtype, "constant")
                plot(xf{i}',"DisplayName", "$\"+fields{1}+"_k="+num2str(vars(i))+"$")
            elseif strcmp(legendtype, "variable1")
                plot(xf{i}',"DisplayName", "$\"+fields{1}+"_k=0.99^k\cdot"+num2str(vars(i))+"$")
            elseif strcmp(legendtype, "variable2")
                plot(xf{i}',"DisplayName", "$\"+fields{1}+"_k=\frac{"+num2str(vars(i))+"}{k}$")
            elseif strcmp(legendtype, "variable3")
                plot(xf{i}',"DisplayName", "$\"+fields{1}+"_k=2^{-length(l_0)}\cdot"+num2str(vars(i))+"$")
            elseif strcmp(legendtype, "nesterov")
                plot(xf{i}',"DisplayName", "$\"+fields{1}+"_k="+num2str(vars(i))+",\;\"+fields{2}+"="+num2str(vars2(i))+"$")
            elseif strcmp(legendtype, "dual_comp")
                if i == 1
                    plot(xf{i}',"DisplayName", "$\"+fields{1}+"_k="+num2str(vars(i))+"$")
                elseif i == 2
                    plot(xf{i}',"DisplayName", "$\"+fields{1}+"_k=0.99^k\cdot"+num2str(vars(i))+"$")
                elseif i == 3
                    plot(xf{i}',"DisplayName", "$\"+fields{1}+"_k=\frac{"+num2str(vars(i))+"}{k}$")
                elseif i == 4
                    plot(xf{i}',"DisplayName", "$\"+fields{1}+"_k=2^{-length(l_0)}\cdot"+num2str(vars(i))+"$")
                else
                    plot(xf{i}',"DisplayName", "$\"+fields{1}+"_k="+num2str(vars(i))+",\;\gamma="+num2str(vars2(i))+"$")
                end
            elseif strcmp(legendtype, "consensus")
                plot(xf{i}',"DisplayName", "$\"+fields{1}+"_k="+num2str(vars(i))+",\;\"+fields{2}+"="+num2str(vars2(i))+"$")
            elseif strcmp(legendtype, "ADMM")
                plot(xf{i}',"DisplayName", "$\"+fields{1}+"="+num2str(vars(i))+"$")
            end
        end
        yscale('log')
        legend()
        xlabel("Iterations")
        ylabel("$\frac{||\textbf{x}_{i,f}-\textbf{x}^*_{f}||_2}{||\textbf{x}^*_{f}||_2}$")
        title("Logarithmic error sequence")
        grid on
        set(gcf, 'Position', get(0, 'Screensize'));
    end

    if savefigure
        set(gcf, 'Position', [1000, 100, 500, 400]);
        set(gcf, "Theme", "light");
        filename = "Images\4_" + filename + ".pdf";
        exportgraphics(gcf,filename)
        close(gcf);
    end
end