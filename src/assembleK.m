function K = assembleK(Ke,nele,nnode)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function to assemble stiffness matrix of a bar
    % IN:   Ke: elemental stiffness
    %       nele: number of elements
    %       nnode: number of nodes
    % OUT:  K: stiffness matrix of a bar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    K = sparse(nnode,nnode);
    
    K = K + sparse(1:nele,1:nele,Ke,nnode,nnode);
    K = K + sparse(1:nele,2:nele+1,-Ke,nnode,nnode);
    K = K + sparse(2:nele+1,1:nele,-Ke,nnode,nnode);
    K = K + sparse(2:nele+1,2:nele+1,Ke,nnode,nnode);


end
