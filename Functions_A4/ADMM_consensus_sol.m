function [traj,xf,u,lambda,fval,exitflag,output,lam] = ADMM_consensus_sol(Planes, rho, plotgen)
    % Function for solving the dual decomposition problem
    tic

    % Initialize parameters
    tol = 1e-3; % Tolerance for x_i(T_final)
    dxf = inf; % Initialize the difference as +inf
    % alpha = 0.1;
    iter = 0; % Iteration counter
    dim = Planes(1).dim;
    Tf = Planes(1).Tf;
    xf = zeros(4*dim.nx,3000);
    u  = zeros(Tf*4*dim.nu,3000);
    opts = optimoptions('quadprog', 'Display', 'off');
    xf_con = zeros(4,3000);
    lambda = zeros(16,1);

    while dxf > tol
        iter = iter + 1;
        for i=1:4
            H = Planes(i).H + rho*Planes(i).A_eq'*Planes(i).A_eq;
            % Avoid Hessian not symmetric warning
            H = (H+H')/2;
            h = Planes(i).h + (0.5*lambda(4*(i-1)+1:4*(i))'*Planes(i).A_eq + rho*(Planes(i).b_eq-xf_con(:,iter))'*Planes(i).A_eq)';
            [u_sol,fval,exitflag,output,lam] = quadprog(H,h,Planes(i).A_u,Planes(i).b_u,[],[],[],[],[],opts);
            u((Tf*2)*(i-1)+1:(Tf*2)*i,iter) = u_sol;
            xf(dim.nx*(i-1)+1:dim.nx*i,iter) = Planes(i).A_eq*u_sol + Planes(i).b_eq;
            xf_con(:,iter+1) = xf_con(:,iter+1) + (Planes(i).A_eq*u_sol + Planes(i).b_eq + lambda(dim.nx*(i-1)+1:dim.nx*i)/rho)/4;
        end
        
        lambda = lambda + rho*(xf(:,iter)-repmat(xf_con(:,iter+1),4,1));

        % Check all relative options for convergence
        dxf = max([norm(xf(dim.nx*0+1:dim.nx*1,iter)-xf(dim.nx*1+1:dim.nx*2,iter)), ...
                   norm(xf(dim.nx*0+1:dim.nx*1,iter)-xf(dim.nx*2+1:dim.nx*3,iter)), ...
                   norm(xf(dim.nx*0+1:dim.nx*1,iter)-xf(dim.nx*3+1:dim.nx*4,iter)), ...
                   norm(xf(dim.nx*1+1:dim.nx*2,iter)-xf(dim.nx*2+1:dim.nx*3,iter)), ...
                   norm(xf(dim.nx*2+1:dim.nx*3,iter)-xf(dim.nx*3+1:dim.nx*4,iter)), ...
                   norm(xf(dim.nx*1+1:dim.nx*2,iter)-xf(dim.nx*3+1:dim.nx*4,iter))]);
        dxf = dxf+1;
        % if abs(dxf1 - dxf) < tol*1e-2 || iter == 3000
        if iter == 1000
            break
        end
    end
    toc

    xf = xf(:,1:iter);

    Tf = Planes(1).Tf;
    
    % Calculate the trajectory for each plane
    for i = 1:4
        trajectory = [Planes(i).x0, reshape(Planes(i).T*Planes(i).x0 + Planes(i).S*u((Tf*2)*(i-1)+1:(Tf*2)*i,iter),4,[])]; 
        eval(sprintf('traj.x%d = trajectory;', i));
    end

    %% Plot
    if plotgen
        figure(2), clf;
        hold on
        plot(traj.x1(1,:),traj.x1(2,:))
        plot(traj.x2(1,:),traj.x2(2,:))
        plot(traj.x3(1,:),traj.x3(2,:))
        plot(traj.x4(1,:),traj.x4(2,:))
        % plot(xf(1),xf(2),'Marker','+', 'MarkerSize',15, 'LineWidth',2)
    end
end