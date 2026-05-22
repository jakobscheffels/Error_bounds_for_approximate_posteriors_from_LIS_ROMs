function u = febar( D, F )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function to evaluate displacement of bar for given load 
    % IN:   D: rigidity of bar
    %       F: load vector as column of matrix
    % OUT:  u: displacement of bar 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('Parameters.mat')
    
    % assemble stiffness matrix
    Ke = D/l;
    
    K = assembleK(Ke,nele,nnode);
    K = applyBoundaryCondition(K,BC_dofs);
    
    % solution for nodal displacements
    u = K\F;

end