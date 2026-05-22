function Ke = elementalStiffnessPlate(coords,E,nu)
    % Plane stress constitutive matrix
    D = E/(1-nu^2)*[1 nu 0;
                    nu 1 0;
                    0 0 (1-nu)/2];

    gp = [-1/sqrt(3)  1/sqrt(3)];
    Ke = zeros(8,8);
    for i=1:2
        for j=1:2
            psi = gp(i);
            eta = gp(j);
            
            J = Jacobian(coords(:,1),coords(:,2),psi,eta);
            
    
            B = BMatrix(J,psi,eta);
    
            Ke = Ke + B'*D*B*det(J);
            
        end
    end
end