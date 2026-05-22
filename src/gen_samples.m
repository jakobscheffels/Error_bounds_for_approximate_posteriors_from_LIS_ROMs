function [S_q, S_u, Phi, Sigma] = gen_samples(n_in)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function to generate n samples and perform SVD of snapshots
    % IN:   n: number of samples (150)
    %       
    % OUT:  S_q: (nelexn) matrix containing realizations of distr. load
    %       S_u: ((nele+1)xn) matrix containing respective displacements
    %       U: (nelexnele) matrix containing left singular vectors of SVD
    %       S: (nelexn) diagonal matrix containing singular values
    %       V:(nxn) matrix containing right singular vectors of SVD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    load Parameters.mat
    
    if nargin ==1
        n=n_in;
        
    else
        n=150;
      
    end

    xGL=0;
    XQ = l/2*(repmat(xGL',1,nele)+repmat(1:2:2*nele,length(xGL),1));
    xq = XQ(:);
    %phi_at_xq = cell2mat(cellfun(@(c) c(xq'),phi,'un',0));
    gamma_prior_q=zeros(nele,nele);
    for i = 1:nele
        for j=1:nele
            gamma_prior_q(i,j)=sigma_q^2*exp(-abs(xq(i)-xq(j))/theta);
        end
    end
    S_q=zeros(nele,n);
    S_u=zeros(nele+1,n);
    
    xi = randn(n,nele);

    S_q = (mu_q.*ones(nele,1)+chol(gamma_prior_q,"lower")*xi');
    F = L_mat*S_q;
    F = applyBoundaryCondition(F,BC_dofs,'Force');

    q = mvnrnd(mu_q.*ones(nele,1),gamma_prior_q,1);
        
    S_u = evaluateU(F);

    [Phi, Sigma, ~] = svd(S_u); 
    
end

function u = evaluateU(F)
    load 'Parameters.mat' 'beam_bool' 'nele' 'D' 'I' 'E'  'BC_dofs' 't' 'zeta' 'k_vector' 'l' 'tunnel'
    if beam_bool
        if tunnel
            EI = pi/64*zeta*E*(D^4-(D-2*t)^4);
            K = assembleBeamK(EI,1,l,nele,nele+1);
            KG = assembleGroundStiffnessK(k_vector);
            K = K+KG;
            K = applyBoundaryCondition(K,BC_dofs);
            u=K\F;
        else
            K = assembleBeamK(E,I,l,nele,nele+1);
            K = applyBoundaryCondition(K,BC_dofs);
            u=K\F;
        end
    else
        u = febar(D,F);
    end
end