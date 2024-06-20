function [xf_avg,xf,traj,u] = general_sol(Planes, simtype, xf_central, avg)
    % Function for solving each problem
    % tic

    % Initialize parameters
    tol = simtype.tolerance; % Tolerance for x_i(T_final)
    max_iter = simtype.iterations;
    dxf = inf; % Initialize the difference as +inf
    dxf_err = dxf;
    iter = 0; % Iteration counter
    dim = Planes(1).dim;
    Tf = Planes(1).Tf;
    xf = zeros(4*dim.nx,max_iter);
    u  = zeros(Tf*4*dim.nu,max_iter);
    opts = optimoptions('quadprog', 'Display', 'off');

    % method specific initialization
    phi = simtype.phi;
    alpha = simtype.alpha;

    % Dual
    lambda2 = zeros(4,1);
    lambda3 = zeros(4,1);
    lambda4 = zeros(4,1);
    lambda_D = [lambda2; lambda3; lambda4];
    psi = zeros(size(lambda_D));

    % Consensus
    xf_con = repmat((Planes(1).x0+Planes(2).x0+Planes(3).x0+Planes(4).x0)/4,4,1);
    lambda_con = zeros(16,1);
    W = [0.75, 0.25, 0,    0;
         0.25, 0.5,  0.25, 0;
         0,    0.25, 0.5,  0.25;
         0,    0,    0.25, 0.75];
    W=kron(W,eye(4));

    % ADMM
    xf_con_ADMM = zeros(4,max_iter+1);
    lambda_ADMM = zeros(16,1);
    rho = simtype.rho;

    while dxf > tol
        iter = iter + 1;

        if strcmp(simtype.method, "dual")
            % Dual decomposition with varying step size methods
            for i=1:4
                if i == 1
                    h = Planes(i).h + (0.5*(lambda_D(1:4,iter)+lambda_D(5:8,iter)+lambda_D(9:12,iter))'*Planes(i).A_eq)';
                else
                    h = Planes(i).h - (0.5*lambda_D(4*(i-2)+1:4*(i-1),iter)'*Planes(i).A_eq)';
                end
                [u_sol,~,~,~,~] = quadprog(Planes(i).H,h,Planes(i).A_u,Planes(i).b_u,[],[],[],[],[],opts);
                u((Tf*2)*(i-1)+1:(Tf*2)*i,iter) = u_sol;
                xf(dim.nx*(i-1)+1:dim.nx*i,iter) = Planes(i).A_eq*u_sol + Planes(i).b_eq;
            end
    
            gamma = 0;
            if strcmp(simtype.update_sequence, "constant")
                % Constant step size
                alpha = simtype.alpha;
            elseif strcmp(simtype.update_sequence, "variable1")
                % Variable step size
                alpha = 0.9999^iter*simtype.alpha0;
            elseif strcmp(simtype.update_sequence, "variable2")
                % Variable step size
                alpha = simtype.alpha/((iter+100)/100);
            elseif strcmp(simtype.update_sequence, "variable3")
                % Variable step size
                alpha = simtype.alpha*2^(-countLeadingZeros(dxf_err));
            elseif strcmp(simtype.update_sequence, "nesterov")
                % Nesterov's accelerated gradient
                alpha = simtype.alpha;
                gamma = simtype.gamma;
            end
            gradf = alpha*[xf(1:dim.nx*1,iter)-xf(dim.nx*(2-1)+1:dim.nx*2,iter);
                           xf(1:dim.nx*1,iter)-xf(dim.nx*(3-1)+1:dim.nx*3,iter);
                           xf(1:dim.nx*1,iter)-xf(dim.nx*(4-1)+1:dim.nx*4,iter)];
            psi = gamma*psi + gradf;
            lambda_D(:,iter+1) = lambda_D(:,iter) + psi;

        elseif strcmp(simtype.method, "consensus")
            % Consensus
            for i=1:4
                h = Planes(i).h + (0.5*lambda_con(4*(i-1)+1:4*(i))'*Planes(i).A_eq)';
                [u_sol,~,~,~,~] = quadprog(Planes(i).H,h,Planes(i).A_u,Planes(i).b_u,[],[],[],[],[],opts);
                u((Tf*2)*(i-1)+1:(Tf*2)*i,iter) = u_sol;
                xf(dim.nx*(i-1)+1:dim.nx*i,iter) = Planes(i).A_eq*u_sol + Planes(i).b_eq;
                xf_con(dim.nx*(i-1)+1:dim.nx*i) = Planes(i).A_eq*u_sol + Planes(i).b_eq;
            end
            xf_con = W^phi*(xf_con);
            lambda_con = lambda_con + alpha*(xf(:,iter)-xf_con);

        elseif strcmp(simtype.method, "ADMM")
            % ADMM
            for i=1:4
                H = Planes(i).H + rho*(Planes(i).A_eq'*Planes(i).A_eq);
                % % Avoid Hessian not symmetric warning
                % H = (H+H')/2;
                h = Planes(i).h + (0.5*lambda_ADMM(4*(i-1)+1:4*(i))'*Planes(i).A_eq + rho*(Planes(i).b_eq-xf_con_ADMM(:,iter))'*Planes(i).A_eq)';
                [u_sol,~,~,~,~] = quadprog(H,h,Planes(i).A_u,Planes(i).b_u,[],[],[],[],[],opts);
                u((Tf*2)*(i-1)+1:(Tf*2)*i,iter) = u_sol;
                xf(dim.nx*(i-1)+1:dim.nx*i,iter) = Planes(i).A_eq*u_sol + Planes(i).b_eq;
                xf_con_ADMM(:,iter+1) = xf_con_ADMM(:,iter+1) + (Planes(i).A_eq*u_sol + Planes(i).b_eq + lambda_ADMM(dim.nx*(i-1)+1:dim.nx*i)/rho)/4;
            end
            lambda_ADMM = lambda_ADMM + rho*(xf(:,iter)-repmat(xf_con_ADMM(:,iter+1),4,1));
        end


        % Check all relative options for convergence
        dxf_err = max([norm(xf(dim.nx*0+1:dim.nx*1,iter)-xf(dim.nx*1+1:dim.nx*2,iter)), ...
                   norm(xf(dim.nx*0+1:dim.nx*1,iter)-xf(dim.nx*2+1:dim.nx*3,iter)), ...
                   norm(xf(dim.nx*0+1:dim.nx*1,iter)-xf(dim.nx*3+1:dim.nx*4,iter)), ...
                   norm(xf(dim.nx*1+1:dim.nx*2,iter)-xf(dim.nx*2+1:dim.nx*3,iter)), ...
                   norm(xf(dim.nx*2+1:dim.nx*3,iter)-xf(dim.nx*3+1:dim.nx*4,iter)), ...
                   norm(xf(dim.nx*1+1:dim.nx*2,iter)-xf(dim.nx*3+1:dim.nx*4,iter))]);
        if strcmp(simtype.endcondition, "convergence")
            dxf = dxf_err;
        end
        if iter == max_iter
            break
        end
    end
    % toc

    xf = xf(:,1:iter);
    Tf = Planes(1).Tf;
    xf_avg = err_norm(xf_central,xf,avg);
    
    % Calculate the trajectory for each plane
    for i = 1:4
        trajectory = [Planes(i).x0, reshape(Planes(i).T*Planes(i).x0 + Planes(i).S*u((Tf*2)*(i-1)+1:(Tf*2)*i,iter),4,[])]; 
        eval(sprintf('traj.x%d = trajectory;', i));
    end
end

function numLeadingZeros = countLeadingZeros(number)
    % Convert the number to a string
    numberStr = num2str(number, '%.15f'); % Ensure sufficient precision
    
    % Find the position of the decimal point
    decimalIndex = find(numberStr == '.', 1);
    
    % Initialize the count of leading zeros
    numLeadingZeros = 0;
    if any(decimalIndex ~= 2) || ~strcmp(numberStr(1),"0")
        return
    end
    
    % Loop through the characters after the decimal point
    for i = decimalIndex+1:length(numberStr)
        if numberStr(i) == '0'
            numLeadingZeros = numLeadingZeros + 1;
        else
            break;
        end
    end
end