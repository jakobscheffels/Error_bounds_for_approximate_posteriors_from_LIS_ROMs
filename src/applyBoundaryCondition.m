function [A]=applyBoundaryCondition(A,dofs,mode)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function to apply boundary condition on matrix A
    % IN:   A: matrix/vector to apply boundary condition to
    %       dofs: fixed degrees of freedom
    %       mode: mode condition is applied to A (Observation, Coordinate,
    %       Index, Force, default)
    %
    % OUT:  A: matrix/vector with boundary condition applied
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('Parameters.mat','beam_bool','fixed')
    counter = 0;
    if size(dofs,1)==0
        return;
    end

    if nargin > 2
        
        if strcmp(mode,'Observation')
            for i=1:size(dofs,2)
                if beam_bool
                    if strcmp(fixed(i),'true')
                        column = 2*dofs(i)-1-counter;
                        A(:,column)=[];
                        A(:,column)=[];
                        counter = counter + 2;
                    else
                        column = 2*dofs(i)-1-counter;
                        A(:,column)=[];
                        counter = counter +1;
                    end
                else
                    column = dofs(i)-(i-1);
                    A(:,column)=[];
                end
            end
        elseif strcmp(mode,'Coordinate')
            x_dofs = A;
            rot_dofs = A;

            for i=1:size(dofs,2)
                column = dofs(i)-(i-1);
                x_dofs(:,column)=[];
            end
            A = x_dofs;
            if beam_bool
                counter = 0;
                for i=1:size(dofs,2)
                    if strcmp(fixed(i),'true')
                        column = dofs(i)-counter;
                        rot_dofs(:,column)=[];
                        counter = counter +1;
                    end
                end
                A = {A,rot_dofs};
            end
        elseif strcmp(mode,'Index')
            if beam_bool
                index_disp = A(1:2:end);
                index_rot = A(2:2:end);
            else
                index_disp = A;
            end

            for i=1:size(dofs,2)
                column = dofs(i)-(i-1);
                index_disp(:,column)=[];
            end
            A = index_disp;
            if beam_bool
                counter = 0;
                for i=1:size(dofs,2)
                    if strcmp(fixed(i),'true')
                        column = dofs(i)-counter;
                        index_rot(:,column)=[];
                        counter = counter +1;
                    end
                end

                counter = 1;
                counter_disp=1;
                counter_rot=1;
                index_disp_it=[];
                index_rot_it=[];
                bool = true;
                while bool
                    if index_disp(counter_disp)<index_rot(counter_rot)
                        index_disp_it(end+1)=counter;
                        counter_disp = counter_disp+1;
                    else
                        index_rot_it(end+1)=counter;
                        counter_rot = counter_rot+1;
                    end
                    counter = counter+1;
                    if counter_disp > size(index_disp,2)
                        index_rot_it = [index_rot_it, counter:(counter+(size(index_rot,2)-counter_rot))];
                        bool = false;
                    elseif counter_rot > size(index_rot,2)
                        index_disp_it = [index_disp_it, counter:(counter+(size(index_disp,2)-counter_disp))];
                        bool = false;
                    end
                end
                
                index_disp = index_disp_it;
                index_rot = index_rot_it;


                A = {index_disp,index_rot};
            end

        elseif strcmp(mode,'Force')
            for i=1:size(dofs,2)
                if beam_bool
                    if strcmp(fixed(i),'true')
                        column = 2*dofs(i)-1-counter;
                        A(column,:)=[];
                        A(column,:)=[];
                        counter = counter +2;
                    else
                        column = 2*dofs(i)-counter;
                        A(column,:)=[];
                        counter = counter +1;
                    end
                else
                    column = dofs(i)-(i-1);
                    A(column,:)=[];
                end
            end
        end
        
    else
        for i=1:size(dofs,2)            
            if beam_bool
                if strcmp(fixed(i),'true')
                    column = 2*dofs(i)-1-counter;
                    A (column,:)=[];
                    A(:,column)=[];
                    A(column,:)=[];
                    A(:,column)=[];
                    counter = counter +2;
                else
                    column = 2*dofs(i)-1-counter;
                    A (column,:)=[];
                    A(:,column)=[];
                    counter = counter +1;
                end
            else
                column = dofs(i)+1-i;
                A (column,:)=[];
                A(:,column)=[];
            end
        end
    end
end