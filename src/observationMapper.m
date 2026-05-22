function C = observationMapper(m_pos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Function to assemble mapping matrices
%%%
%%%     Input:  m_pos      Position indices of observable nodes
%%%
%%%     Output: C          Operator to extract observable nodes from full
%%%                         state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load Parameters.mat 'beam_bool' 'nnode'
    
    m = length(m_pos);

    if beam_bool
        C = zeros(m,2*nnode);
        j=1;
        for i=1:m
            C(j,2*(m_pos(i)-1)+1)=1;
            j=j+1;
        end
    else
        C = zeros(m,nnode);
        j=1;
        for i=1:m
            C(j,m_pos(i))=1;
            j=j+1;
        end        
    end
end