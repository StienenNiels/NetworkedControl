function trajectories_plot(xf_central,traj_compare,filename)
    if nargin < 3
        savefigure = 0;
        plotstate = 1;
    else
        plotstate = 1;
        savefigure = 1;
    end

    if plotstate
        figure(99),clf
        tiledlayout(4,1)
        for i = 1:4
            nexttile
            if i == 1
                title("Aircraft trajectories")
                ylabel("$x\; [m]$")
            elseif i == 2
                ylabel("$y\; [m]$")
            elseif i == 3
                ylabel("$\dot{x}\; [m/s]$")
            elseif i == 4
                ylabel("$\dot{y}\; [m/s]$")
            end
            hold on
            plot(0:1:5, traj_compare.x1(i,:))
            plot(0:1:5, traj_compare.x2(i,:))
            plot(0:1:5, traj_compare.x3(i,:))
            plot(0:1:5, traj_compare.x4(i,:))
            yline(xf_central(i), "LineStyle","-.","LineWidth",1.5)
            xlim([0,7])
            grid on
            legend("Plane 1", "Plane 2", "Plane 3", "Plane 4", "Centralized", Location="east")
        end
        xlabel("Time")
    end

    if savefigure
        set(gcf, 'Position', [1000, 100, 500, 600]);
        set(gcf, "Theme", "light");
        filename = "Images\4_" + filename + ".pdf";
        exportgraphics(gcf,filename)
        close(gcf);
    end
end