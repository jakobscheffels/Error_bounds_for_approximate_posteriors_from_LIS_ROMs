function [posCov, G, df] = solveReducedModel(V,W)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function to solve for posterior covariance using model reduction
    % IN:   V: reconstruction basis 
    %       W: projection basis
    % OUT:  posCov: posterior covariance approximation
    %       G: approximation of forward operator using reduced matrix
    %       df: foerstner distance to true posterior covariance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load LIP_Setup.mat K C
    
    % Build reduced stiffness matrix
    K_red = W'*K*V;
    
    % Construct approximation of forward operator
    G = C*V*(K_red\W');
    
    % Solve for posterior covariance approximation and Foerstner distance
    [posCov, df] = posteriorCovarianceReduced(G);
end