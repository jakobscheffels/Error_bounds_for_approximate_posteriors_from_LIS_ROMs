function Parameters(inputName, inputVal)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     Input: inputName is cell containing variable names to be modified
    %%%            inputVal is list containing values of variables
    %%%
    %%%     Output: Stores modified parameters in Parameter.mat to be loaded
    %%%
    %%%     Valid variable names (default value):
    %%%       mu_q: mean of distributed load (1)
    %%%       sigma_q: standard deviation of distributed load (0.2)
    %%%       L: length of bar element (2)
    %%%       theta: correlation length (2)
    %%%       D: rigidity of bar (4e8)
    %%%       E: Youngs modulus of beam (10000)
    %%%       nele: number of finite elements used (100)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % scale of fluctuation of q (=2*correlation length)
    theta = 2;
    % number of elements
    nele = 100;
    nnode = nele+1;
    % length of domain
    L = 2;

    D = 4e8;
    E = 10000;
    I = 10000;
    % DOF with Dirichlet B.C.
    BC_dofs = 1;
    
    % Boolean to determine to use beam elements
    beam_bool = false;
    % Dirichlet B.C.
    fixed = {'true'};

    % distributed load
    mu_q = 4e6;
    sigma_q = 0.3*mu_q;
    
    % Boolean to use tunnel segments
    tunnel = false;
    zeta = 1/7;
    t = 0.35;
    % set Parameters
    if nargin > 0
        for i=1:size(inputName,2)
            input = inputVal(i);
            if strcmp(inputName(i),'theta')
                theta=input{1};      
            elseif strcmp(inputName(i),'nele')
                nele=input{1};
                nnode=nele+1;
            elseif strcmp(inputName(i),'L')
                L=input{1};
            elseif strcmp(inputName(i),'D')
                D=input{1};
            elseif strcmp(inputName(i),'E')
                E=input{1};
            elseif strcmp(inputName(i),'I')
                I=input{1};
            elseif strcmp(inputName(i),'mu_q')
                mu_q=input{1};
            elseif strcmp(inputName(i),'sigma_q')
                sigma_q=input{1};
            elseif strcmp(inputName(i),'beam_bool')
                beam_bool=input{1};
            elseif strcmp(inputName(i),'fixed')
                fixed=input{1};
            elseif strcmp(inputName(i),'BC_dofs')
                BC_dofs=sort(input{1});
            elseif strcmp(inputName(i),'tunnel')
                tunnel=input{1};
            elseif strcmp(inputName(i),'k_vector')
                k_vector=input{1};
            elseif strcmp(inputName(i),'zeta')
                zeta=input{1};
            elseif strcmp(inputName(i),'t')
                t=input{1};
            else
                inputStr = inputName(i);
                inputStr = inputStr{:};
                text = append('Invalid variable name: ',inputStr,'! Using default parameter instead!');
                warning(text)
            end
        end
    end
    if tunnel
        nele = 800;
        nnode = nele+1;
        L = 200;
        E = 35e9;
        D = 6.2;
        nele = 800;
        
        BC_dofs=[];
        theta = L/2;
        mu_q = 300e5;
        delta_q = 1;
        sigma_q = mu_q*delta_q;
    end

    k1=33000e3;
    k2=5000e3;

    k_vector = zeros(nele+1,1);
    k_vector(1:nele/2)=k1.*ones(nele/2,1);
    k_vector(nele/2+1:end)=k2.*ones(nele/2+1,1);
    
    l=L/nele;
    x_dofs=0:l:L;
    
    if beam_bool
        index_disp = 1:2*nnode;
    else
        index_disp = 1:nnode;
    end

    L_mat = loadMapper(beam_bool,nnode,nele,l);

    % Store parameters in 'Parameters.mat'
    save('Parameters','BC_dofs','D', ...
        'E','I','l','L','L_mat','mu_q','nele','nnode','sigma_q','theta','tunnel',...
        'k_vector','x_dofs','beam_bool','fixed','index_disp',...
        'zeta','t')
end