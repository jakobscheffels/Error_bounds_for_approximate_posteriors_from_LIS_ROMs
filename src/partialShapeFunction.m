function grad = partialShapeFunction(i,psi,eta)
    
    grad = zeros(2,1);
    if i==1
        grad(1) = 0.25*(1+eta);
        grad(2) = 0.25*(1+psi);
    elseif i==2
        grad(1) = -0.25*(1+eta);
        grad(2) = 0.25*(1-psi);
    elseif  i==3
        grad(1) = -0.25*(1-eta);
        grad(2) = -0.25*(1-psi);
    else
        grad(1) = 0.25*(1-eta);
        grad(2) = -0.25*(1+psi);
    end        
end