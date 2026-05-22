function [delta,V,W]=calculateLISBasisIP(input)
    if nargin > 0
        gamma = input;
    else
        gamma = 1e7;
    end
    load LIP_Setup.mat C K S_pr gamma_obs
    
    W_1 = null(C,'r');
    L = C'/(C*C');
    gamma_obs_eta = gamma_obs(1,1).*eye(size(W_1,2));
    gamma_obs_tilde = L*gamma_obs*L'+gamma^2.*W_1*gamma_obs_eta*W_1';
    S_obs_tilde = sqrt(gamma_obs_tilde);

    R = (K\eye(size(K)))'/S_obs_tilde';
    [U,delta,Z]=svd(R'*S_pr);
    delta = diag(delta);
    tol = max(size(R'*S_pr))*eps(max(delta));
    r = sum(delta>tol);
    U = U(:,1:r);
    delta=delta(1:r);
    Z = Z(:,1:r);
    %C_inv = C'/(C*C');
    
    %V = S_pr*(Z.*(1./sqrt(delta))');
    %W = R*(U.*(1./sqrt(delta))');
    V = S_obs_tilde*U;
    W = R*U.*(1./delta)';

end