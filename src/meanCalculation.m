function mu = meanCalculation(G,y)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate posterior mean
    % IN:   G: forward operator of inverse problem
    %       y: measurement sample
    % OUT:  mu: posterior mean
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load LIP_Setup.mat 'mu_f' 'gamma_prior_f' 'gamma_obs'
    
    % Analyitcal posterior mean for given G and measurement y
    mu = mu_f+gamma_prior_f*G'*((G*gamma_prior_f*G'+gamma_obs)\(y-G*mu_f));
end