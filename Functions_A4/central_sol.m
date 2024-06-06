function [traj,xf,u,fval,exitflag,output,lambda] = central_sol(Planes, plotgen)
    tic
    H=[];h=[];A_eq=[];b_eq=[];A_u=[];b_u=[];
    for i = 1:4
        [H, h, A_eq, b_eq, A_u, b_u] = central_gen(Planes(i), H, h, A_eq, b_eq, A_u, b_u);
    end
    Tfinal = Planes(1).Tf;

    % Solve the central optimization problem
    opts = optimoptions('quadprog', 'Display', 'off');
    [u,fval,exitflag,output,lambda] = quadprog(H,h,A_u,b_u,A_eq,-b_eq,[],[],[],opts);
    
    % Calculate the trajectory for each plane
    for i = 1:4
        trajectory = [Planes(i).x0, reshape(Planes(i).T*Planes(i).x0 + Planes(i).S*u((Tfinal*2)*(i-1)+1:(Tfinal*2)*i),4,[])]; 
        eval(sprintf('traj.x%d = trajectory;', i));
    end
    % Final position where planes are converged
    xf = traj.x1(:,end);

    % Spruce it up a bit bro, you can do better
    if plotgen
        figure(1),clf
        hold on
        plot(traj.x1(1,:),traj.x1(2,:))
        plot(traj.x2(1,:),traj.x2(2,:))
        plot(traj.x3(1,:),traj.x3(2,:))
        plot(traj.x4(1,:),traj.x4(2,:))
        plot(xf(1),xf(2),'Marker','+', 'MarkerSize',15, 'LineWidth',2)
    end
    toc
end

% Combine all matrices into one large problem to solve with quadprog
function [H, h, A_eq, b_eq, A_u, b_u] = central_gen(plane, H, h, A_eq, b_eq, A_u, b_u)
    dim = plane.dim;
    H = blkdiag(H,plane.H);
    h = [h,plane.h];
    if plane.plane == 1
        A_eq = [repmat(plane.A_eq,[3, 1]),zeros(3*dim.nx,3*dim.nu*dim.N)];
        b_eq = repmat(plane.b_eq,[3, 1]);
    else
        id = plane.plane;
        A_eq((id-2)*dim.nx+1:(id-1)*dim.nx,(id-1)*dim.nu*dim.N+1:(id)*dim.nu*dim.N) = -plane.A_eq;
        b_eq((id-2)*dim.nx+1:(id-1)*dim.nx,:) = b_eq((id-2)*dim.nx+1:(id-1)*dim.nx,:)-plane.b_eq;
    end
    A_u  = blkdiag(A_u,plane.A_u);
    b_u  = [b_u; plane.b_u];
end