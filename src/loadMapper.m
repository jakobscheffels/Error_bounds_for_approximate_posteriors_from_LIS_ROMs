function L_mat = loadMapper(beam_bool,nnode,nele,l)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Function to assemble mapping matrices to map distributed load to
%%%     nodal force vector
%%%
%%%     Input:  beam_bool   True if structure is beam
%%%             nnode       Number of nodes
%%%             nele        Number of elements
%%%             l           Element length
%%%
%%%     Output: L_mat       Operator to map from distributed load to nodal
%%%                         force vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    if beam_bool
        L_mat = zeros(2*nnode,nele);
        for i=1:nele
            L_mat(2*(i-1)+1:2*(i-1)+4,i)=[l/2 l^2/12 l/2 -l^2/12]';
        end
        
    else
        
        L_mat = sparse(nnode,nele);
        L_mat = L_mat + sparse(1:nele,1:nele,l/2,nnode,nele);
        L_mat = L_mat + sparse(2:nele+1,1:nele,l/2,nnode,nele);
    end
end