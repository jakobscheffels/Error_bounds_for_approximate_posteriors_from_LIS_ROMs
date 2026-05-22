function KG = assembleGroundStiffnessK(k_vector)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function to assemble stiffness due to ground support
    % IN:   k_vector: vector containing stiffness of ground
    %       
    % OUT:  KG: stiffness matrix due to ground support
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load 'Parameters.mat'
    
    Dl = D*l;
    KG = sparse(2*nnode,2*nnode);
    for i=1:nele
        k1 = k_vector(i);
        k2 = k_vector(i+1);
        Ke = Dl.*[13/35*k1+3/35*(k2-k1), 11/210*k1*l+1/60*(k2-k1)*l, 9/70*k1+9/140*(k2-k1), -13/420*k1*l-1/70*(k2-k1)*l;...
            11/210*k1*l+1/60*(k2-k1)*l, 1/105*k1*l^2+1/280*(k2-k1)*l^2, 13/420*k1*l+1/60*(k2-k1)*l, -1/140*k1*l^2-1/280*(k2-k1)*l^2;...
            9/70*k1+9/140*(k2-k1), 13/420*k1*l+1/60*(k2-k1)*l, 13/35*k1+2/7*(k2-k1), -11/210*k1*l-1/28*(k2-k1)*l;...
            -13/420*k1*l-1/70*(k2-k1)*l, -1/140*k1*l^2-1/280*(k2-k1)*l^2, -11/210*k1*l-1/28*(k2-k1)*l, 1/105*k1*l^2+1/168*(k2-k1)*l^2];
        KG(2*(i-1)+1:(2*(i-1)+4),2*(i-1)+1:(2*(i-1)+4))=KG(2*(i-1)+1:(2*(i-1)+4),2*(i-1)+1:(2*(i-1)+4))+Ke;
    
    end
end