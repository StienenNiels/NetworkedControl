function [traj,xf,u,lambda,fval,exitflag,output,lam] = consensus_sol(Planes, alpha, phi, plotgen)
    % Function for solving the dual decomposition problem
    tic

    % Initialize parameters
    tol = 1e-3; % Tolerance for x_i(T_final)
    dxf = inf; % Initialize the difference as +inf
    % alpha = 0.1;
    iter = 0; % Iteration counter
    dim = Planes(1).dim;
    Tf = Planes(1).Tf;
    xf = zeros(4*dim.nx,2700);
    u  = zeros(Tf*4*dim.nu,2700);
    opts = optimoptions('quadprog', 'Display', 'off');
    xf_con = repmat((Planes(1).x0+Planes(2).x0+Planes(3).x0+Planes(4).x0)/4,4,1);
    lambda = zeros(16,1);
    W = [0.75, 0.25, 0,    0;
         0.25, 0.5,  0.25, 0;
         0,    0.25, 0.5,  0.25;
         0,    0,    0.25, 0.75];
    W=kron(W,eye(4));

    while dxf > tol
        iter = iter + 1;
        for i=1:4
            h = Planes(i).h + (0.5*lambda(4*(i-1)+1:4*(i))'*Planes(i).A_eq)';
            [u_sol,fval,exitflag,output,lam] = quadprog(Planes(i).H,h,Planes(i).A_u,Planes(i).b_u,[],[],[],[],[],opts);
            u((Tf*2)*(i-1)+1:(Tf*2)*i,iter) = u_sol;
            xf_con(dim.nx*(i-1)+1:dim.nx*i) = Planes(i).A_eq*u_sol + Planes(i).b_eq;
            xf(dim.nx*(i-1)+1:dim.nx*i,iter) = xf_con(dim.nx*(i-1)+1:dim.nx*i);
        end
        
        xf_con = W^phi*(xf_con);
        lambda = lambda + alpha*(xf(:,iter)-xf_con);

        % Check all relative options for convergence
        dxf = max([norm(xf(dim.nx*0+1:dim.nx*1,iter)-xf(dim.nx*1+1:dim.nx*2,iter)), ...
                   norm(xf(dim.nx*0+1:dim.nx*1,iter)-xf(dim.nx*2+1:dim.nx*3,iter)), ...
                   norm(xf(dim.nx*0+1:dim.nx*1,iter)-xf(dim.nx*3+1:dim.nx*4,iter)), ...
                   norm(xf(dim.nx*1+1:dim.nx*2,iter)-xf(dim.nx*2+1:dim.nx*3,iter)), ...
                   norm(xf(dim.nx*2+1:dim.nx*3,iter)-xf(dim.nx*3+1:dim.nx*4,iter)), ...
                   norm(xf(dim.nx*1+1:dim.nx*2,iter)-xf(dim.nx*3+1:dim.nx*4,iter))]);
        dxf = dxf+1;
        % if abs(dxf1 - dxf) < tol*1e-2 || iter == 3000
        if iter == 100
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