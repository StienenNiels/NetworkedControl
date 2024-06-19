function [traj,xf,u,lambda,fval,exitflag,output,lam] = dual_sol(Planes, alpha0, update_seq, plotgen)
    % Function for solving the dual decomposition problem
    tic

    % Initialize parameters
    tol = 1e-3; % Tolerance for x_i(T_final)
    dxf = inf; % Initialize the difference as +inf
    iter = 0; % Iteration counter
    dim = Planes(1).dim;
    Tf = Planes(1).Tf;
    xf = zeros(4*dim.nx,2700);
    u  = zeros(Tf*4*dim.nu,2700);
    opts = optimoptions('quadprog', 'Display', 'off');
    lambda2 = zeros(4,1);
    lambda3 = zeros(4,1);
    lambda4 = zeros(4,1);
    lambda = [lambda2; lambda3; lambda4];
    phi = zeros(size(lambda));
    dxf1 = 0;

    while dxf > tol
        iter = iter + 1;
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

        if update_seq == 0
            % Constant step size
            alpha = alpha0;
            gamma = 0;
        elseif update_seq == 1
            % Variable step size
            alpha = 0.9999^iter*alpha0;
            gamma = 0;
        elseif update_seq == 2
            % Nesterov's accelerated gradient
            alpha = alpha0;
            gamma = 0.8;
        end

        gradf = alpha*[xf(1:dim.nx*1,iter)-xf(dim.nx*(2-1)+1:dim.nx*2,iter);
                       xf(1:dim.nx*1,iter)-xf(dim.nx*(3-1)+1:dim.nx*3,iter);
                       xf(1:dim.nx*1,iter)-xf(dim.nx*(4-1)+1:dim.nx*4,iter)];

        phi = gamma*phi + gradf;
        lambda(:,iter+1) = lambda(:,iter) + phi;


        % Check all relative options for convergence
        dxf = max([norm(xf(dim.nx*0+1:dim.nx*1,iter)-xf(dim.nx*1+1:dim.nx*2,iter)), ...
                   norm(xf(dim.nx*0+1:dim.nx*1,iter)-xf(dim.nx*2+1:dim.nx*3,iter)), ...
                   norm(xf(dim.nx*0+1:dim.nx*1,iter)-xf(dim.nx*3+1:dim.nx*4,iter)), ...
                   norm(xf(dim.nx*1+1:dim.nx*2,iter)-xf(dim.nx*2+1:dim.nx*3,iter)), ...
                   norm(xf(dim.nx*2+1:dim.nx*3,iter)-xf(dim.nx*3+1:dim.nx*4,iter)), ...
                   norm(xf(dim.nx*1+1:dim.nx*2,iter)-xf(dim.nx*3+1:dim.nx*4,iter))]);
        dxf = dxf+1;
        % if abs(dxf1 - dxf) < tol*1e-2 || iter == 3000
        if iter == 3000
            break
        end
        % abs(dxf1 - dxf)
        dxf1 = dxf;
    end
    toc

    xf = xf(:,1:iter);
    Tf = Planes(1).Tf;
    
    % Calculate the trajectory for each plane
    for i = 1:4
        trajectory = [Planes(i).x0, reshape(Planes(i).T*Planes(i).x0 + Planes(i).S*u((Tf*2)*(i-1)+1:(Tf*2)*i,iter),4,[])]; 
        eval(sprintf('traj.x%d = trajectory;', i));
    end
end