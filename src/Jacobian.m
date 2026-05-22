function J = Jacobian(x,y,eta,psi)
    J = zeros(2,2);
    
    for i=1:4
        grad = partialShapeFunction(i,psi,eta);
        J(1,1) = J(1,1)+grad(1)*x(i);
        J(1,2) = J(1,2)+grad(1)*y(i);
        J(2,1) = J(2,1)+grad(2)*x(i);
        J(2,2) = J(2,2)+grad(2)*y(i);
    end    
end