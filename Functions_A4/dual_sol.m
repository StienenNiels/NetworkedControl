function [traj,xf,u,lambda,fval,exitflag,output,lam] = dual_sol(Planes, plotgen)
    % Function for solving the dual decomposition problem
    tic

    % Initialize parameters
    tol = 1e-3; % Tolerance for x_i(T_final)
    dxf = inf; % Initialize the difference as +inf
    alpha = 0.1;
    iter = 1; % Iteration counter
    dim = Planes(1).dim;
    Tf = Planes(1).Tf;
    xf = zeros(4*dim.nx,2700);
    u  = zeros(Tf*4*dim.nu,2700);
    opts = optimoptions('quadprog', 'Display', 'off');
    lambda2 = zeros(4,1);
    lambda3 = zeros(4,1);
    lambda4 = zeros(4,1);
    lambda = [lambda2; lambda3; lambda4];

    while dxf > tol
        for i=1:4
            if i == 1
                h = Planes(i).h + (0.5*(lambda(1:4,iter)+lambda(5:8,iter)+lambda(9:12,iter))'*Planes(i).A_eq)';
            else
                h = Planes(i).h - (0.5*lambda(4*(i-2)+1:4*(i-1),iter)'*Planes(i).A_eq)';
            end
            [u_sol,fval,exitflag,output,lam] = quadprog(Planes(i).H,h,Planes(i).A_u,Planes(i).b_u,[],[],[],[],[],opts);
            u((Tf*2)*(i-1)+1:(Tf*2)*i,iter) = u_sol;
            xf(dim.nx*(i-1)+1:dim.nx*i,iter) = Planes(i).A_eq*u_sol + Planes(i).b_eq;
        end
        lambda(:,iter+1) = lambda(:,iter) + alpha*[xf(1:dim.nx*1,iter)-xf(dim.nx*(2-1)+1:dim.nx*2,iter);
                                                   xf(1:dim.nx*1,iter)-xf(dim.nx*(3-1)+1:dim.nx*3,iter);
                                                   xf(1:dim.nx*1,iter)-xf(dim.nx*(4-1)+1:dim.nx*4,iter)];
        
        dxf = max([norm(xf(1:dim.nx*1,iter)-xf(dim.nx*1+1:dim.nx*2,iter)), ...
                   norm(xf(1:dim.nx*1,iter)-xf(dim.nx*2+1:dim.nx*3,iter)), ...
                   norm(xf(1:dim.nx*1,iter)-xf(dim.nx*3+1:dim.nx*4,iter))]);
        
        iter = iter + 1;
    end
    toc

    Tf = Planes(1).Tf;
    
    % Calculate the trajectory for each plane
    for i = 1:4
        trajectory = [Planes(i).x0, reshape(Planes(i).T*Planes(i).x0 + Planes(i).S*u((Tf*2)*(i-1)+1:(Tf*2)*i,end),4,[])]; 
        eval(sprintf('traj.x%d = trajectory;', i));
    end
    % Final position where planes are converged
    xf = traj.x1(:,end);

    %% Plot
    % Spruce it up a bit bro, you can do better
    if plotgen
        figure(2),clf
        hold on
        plot(traj.x1(1,:),traj.x1(2,:))
        plot(traj.x2(1,:),traj.x2(2,:))
        plot(traj.x3(1,:),traj.x3(2,:))
        plot(traj.x4(1,:),traj.x4(2,:))
        plot(xf(1),xf(2),'Marker','+', 'MarkerSize',15, 'LineWidth',2)
    end
end