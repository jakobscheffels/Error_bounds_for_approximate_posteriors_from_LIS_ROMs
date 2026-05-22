function [delta,V,W]=calculateLISBasis()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function to calculate LIS basis
    % OUT:  delta: singular values of basis
    %       V: reconstruction basis
    %       W: projection basis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load Parameters.mat BC_dofs
    load LIP_Setup.mat

    R = G'.*(1./sqrt(diag(gamma_obs)))';
    [U,delta,Z]=svd(R'*S_pr);
    delta = diag(delta);
    tol = max(size(R'*S_pr))*eps(max(delta));
    r = sum(delta>tol);
    U = U(:,1:r);
    delta=delta(1:r);
    Z = Z(:,1:r);
    V = S_pr*Z;
    W = R*U.*(1./delta)';
    save LIS_Basis W V 
end