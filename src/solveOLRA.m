function [posCov,G, df] = solveOLRA(V,W)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function to solve for posterior covariance using OLR
    % IN:   V: reconstruction basis
    %       W: projection basis
    % OUT:  posCov: posterior covariance approximation
    %       G: approximation of forward operator
    %       df: foerstner distance to true posterior covariance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load LIP_Setup.mat G

    % Construct optimal low-rank Projector Pr
    Pr = zeros(size(V,1),size(W,1));
    for i=1:size(V,2)
        Pr=Pr+V(:,i)*W(:,i)';
    end
    % Define OLR forward operator approximation
    G = G*Pr;
    
    % Solve for posterior covariance approximation and Foerstner distance
    [posCov,df] = posteriorCovarianceReduced(G);
end