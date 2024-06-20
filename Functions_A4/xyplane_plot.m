function xyplane_plot(xf_central,traj_compare,filename)
    if nargin < 3
        savefigure = 0;
        plotstate = 1;
    else
        plotstate = 1;
        savefigure = 1;
    end

    if plotstate
        figure(99),clf
        tiledlayout(1,1)
        nexttile
        hold on
        plot(traj_compare.x1(1,:),traj_compare.x1(2,:))
        plot(traj_compare.x2(1,:),traj_compare.x2(2,:))
        plot(traj_compare.x3(1,:),traj_compare.x3(2,:))
        plot(traj_compare.x4(1,:),traj_compare.x4(2,:))
        plot(xf_central(1),xf_central(2),'Marker','+','Color','k' , 'MarkerSize',15, 'LineWidth',2)
        % ylim([1e-10, 1e1])
        % % xlim([0,1500])
        legend("Plane 1", "Plane 2", "Plane 3", "Plane 4", "Centralized solution", Location="east")
        xlabel("$x\; [m]$")
        ylabel("$y\; [m]$")
        title("XY plane of aircraft coordination")
        axis("equal")
        grid on
    end

    if savefigure
        % set(gcf, 'Position', [1000, 100, 500, 400]);
        set(gcf, "Theme", "light");
        filename = "Images\4_" + filename + ".pdf";
        exportgraphics(gcf,filename)
        close(gcf);
    end
end