function K = assembleBeamK(E,I,l,nele,nnode)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function to assemble stiffness matrix of beam
    % IN:   E: Youngs modulus
    %       I: bending moment of inertia
    %       l: element length
    %       nele: number of elements
    %       nnode: number of nodes
    %
    % OUT:  K: stiffness matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EI = E*I/l^3;
    Ke =EI.*[12 6*l -12 6*l; 6*l 4*l^2 -6*l 2*l^2;-12 -6*l 12 -6*l; 6*l 2*l^2 -6*l 4*l^2];
    K = sparse(2*nnode,2*nnode);
    
    for i = 1:nele
        K(2*(i-1)+1:(2*(i-1)+4),2*(i-1)+1:(2*(i-1)+4))=K(2*(i-1)+1:(2*(i-1)+4),2*(i-1)+1:(2*(i-1)+4))+Ke;
    end
end