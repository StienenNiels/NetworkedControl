function simtype = opt_sim(method, update_sequence, alpha, gamma, rho, phi, endcondition, iterations, tolerance)
    if strcmp(method,"ADMM") || strcmp(method,"dual") || strcmp(method,"consensus")
        simtype.method = method;
    else
        error("Invalid method chosen");
    end
    if strcmp(endcondition,"iteration") || strcmp(endcondition,"convergence")
        simtype.endcondition = endcondition;
    else
        error("Invalid endcondition chosen");
    end
    if strcmp(update_sequence,"constant") || strcmp(update_sequence,"variable1") || ...
            strcmp(update_sequence,"variable2") || strcmp(update_sequence,"variable3") || ...
            strcmp(update_sequence,"nesterov")
        simtype.update_sequence = update_sequence;
    else
        error("Invalid update sequence chosen");
    end
    simtype.alpha = alpha;
    simtype.alpha0 = alpha;
    simtype.phi = phi;
    simtype.rho = rho;
    simtype.gamma = gamma;
    simtype.iterations = iterations;
    simtype.tolerance = tolerance;
end