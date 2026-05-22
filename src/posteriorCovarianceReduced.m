function [posCov, df] = posteriorCovarianceReduced(G)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function to evaluate posterior covariance and foersnter distance
    % IN:   G: forward operator of inverse problem
    %       
    % OUT:  posCov: posterior covariance for given inverse problem
    %       df: foerstner distance to true posterior covariance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load LIP_Setup.mat gamma_obs gamma_prior_f_inv singular_prior gamma_prior_f
    load LIS_Basis.mat V W

    prior_hat = gamma_prior_f;
    posCov = prior_hat-prior_hat*G'*((G*prior_hat*G'+gamma_obs)\G)*prior_hat;

    df = foerstnerDistance(posCov);
end

function df = foerstnerDistance(gamma1)
    load LIP_Setup.mat gamma_pos
    load LIS_Basis.mat W
  
    % Calculate cholesky factor of projected covariances
    L = chol(W'*gamma1*W,'lower');
    R = chol(W'*gamma_pos*W,'lower');
    %L = chol(gamma1,'lower');
    %R = chol(gamma_pos,'lower');
    % Solve for generalized eigenvalues
    df=gen_eigenvalue(L,R);
    df(abs(df)<eps)=[];
    % Evaluate squared Foerstner distance
    df = dot(log(df),log(df));
end

function delta = gen_eigenvalue(L,R)
    % Solve for generalized eigenvalues using SVD
    tmp = L'/R';
    [~,S,~]=svd(tmp);
    delta = diag(S).^2;
end