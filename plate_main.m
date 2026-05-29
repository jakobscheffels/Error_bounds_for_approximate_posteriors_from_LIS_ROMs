clear all
addpath src\

rng(42)
load plate_no_hole.mat

nodes = size(coordinatesFEM,1);
num_ele = size(elementsFEM,1);
ndof = 2*nodes;

%% determine dofs of boundary condition
Y_bottom = min(coordinatesFEM(:,2));
tol = 1e-8;
plate_width = max(coordinatesFEM(:,1))-min(coordinatesFEM(:,1));

bottom_nodes = find(abs(coordinatesFEM(:,2) - Y_bottom) < tol);
dofs_bottom = zeros(2*size(bottom_nodes,1),1);
for i=1:size(bottom_nodes,1)
    dofs_bottom(2*i-1) = 2.*bottom_nodes(i)-1;
    dofs_bottom(2*i) = 2.*bottom_nodes(i);
end

Y_top = max(coordinatesFEM(:,2));
top_nodes = find(abs(coordinatesFEM(:,2)-Y_top)<tol);
dofs_top = zeros(2*size(top_nodes,1),1);
dofs_topY = zeros(size(top_nodes,1),1);
for i=1:size(top_nodes,1)
    dofs_top(2*i-1)=2.*top_nodes(i)-1;
    dofs_top(2*i)=2.*top_nodes(i);
    dofs_topY(i)=2.*top_nodes(i);
end
coords_top = coordinatesFEM(top_nodes,:);


%% Loading 
l_c = 0.3*plate_width;
mu=1e6;
delta=0.3;
mu_q = zeros(length(top_nodes),1);

nSegments = length(top_nodes)-1;
len = zeros(nSegments,1);
x_mid = zeros(nSegments,1);
for i=1:nSegments
    x_left = coordinatesFEM(top_nodes(i),1);
    x_right = coordinatesFEM(top_nodes(i+1),1);
    x_mid(i) = (x_right+x_left)/2;
    len(i)=x_right-x_left;
end
sigma = delta*mu;
B = zeros(nSegments,nSegments);
for i=1:nSegments
    for j=1:nSegments
        B(i,j)=sigma^2*exp(-abs(x_mid(i)-x_mid(j))/l_c);
    end
    mu_q(i:i+1)=mu_q(i:i+1)+len(i)/2.*[mu;mu];
end
L = chol(B,"lower");
z=randn(nSegments,1);
q_realization = mu+L*z;

N_rep = 1;
q_mean = mu+L*randn(nSegments,N_rep);

F = zeros(ndof,1);
mu_f = zeros(ndof,1);

for i=1:nSegments
    node_left = top_nodes(i);
    node_right = top_nodes(i+1);
    L_seg = coordinatesFEM(node_right,1)-coordinatesFEM(node_left,1);
    f_local = q_realization(i)*L_seg/2 * [1;1];

    F(dofs_top(2*i))=F(dofs_top(2*i))+f_local(1);
    F(dofs_top(2*(i+1)))=F(dofs_top(2*(i+1)))+f_local(2);

   
    f_local = q_mean(i,:).*L_seg/2.*[1;1];
end
mu_f(dofs_topY)=mu_q;



%% IP Parameters

gamma_obs = 2e-4;

m=6;

% measurement locations at top 
m_pos = randi([1,17],1,m);
m_pos = unique(m_pos);

while ~(size(unique(m_pos),2)==m)
    m_pos=[m_pos, randi([1,17],1,m-size(m_pos,2))];
    m_pos = unique(m_pos);
end

m_pos=sort(m_pos);
m_dofs = dofs_topY(m_pos);
m_dofs = sort(m_dofs);

% observation covariance
gamma_obs = gamma_obs.^2*eye(m);
S_obs = sqrt(gamma_obs);

% prior covariance of discretized load
gamma_prior_q=zeros(nSegments,nSegments);
x_coords_top = coords_top(:,1);
for i = 1:nSegments
    for j=1:nSegments
        gamma_prior_q(i,j)=sigma^2*exp(-abs((x_coords_top(i)+x_coords_top(i+1))/2-((x_coords_top(j)+x_coords_top(j+1))/2))/l_c);
    end
end
S_q = chol(gamma_prior_q,"lower");

% prior covariance of nodal force vector
S_pr = zeros(ndof,size(S_q,1));
S_pr (dofs_topY(1:end-1),:)=L_seg/2.*S_q;
S_pr(dofs_topY(2:end),:) =S_pr(dofs_topY(2:end),:)+L_seg/2.*S_q;


gamma_prior_f = zeros(ndof,ndof);
gamma_prior_f(dofs_topY(1:end-1),dofs_topY(1:end-1))=L_seg/2.*gamma_prior_q;
gamma_prior_f(dofs_topY(2:end),dofs_topY(2:end))=gamma_prior_f(dofs_topY(2:end),dofs_topY(2:end))+L_seg/2.*gamma_prior_q;
gamma_prior_full = gamma_prior_f;
gamma_prior_f(dofs_bottom,:)=[];
gamma_prior_f(:,dofs_bottom)=[];
S_pr (dofs_bottom,:)=[];

C = zeros(m,ndof);
C (:,m_dofs)=eye(m);
C (:,dofs_bottom)=[];

%% Elemental degree of freedom table
elementDofs = zeros(num_ele,8);
for i = 1:num_ele
    idx = elementsFEM(i,:);
    for j=1:4
        elementDofs(i,2*j-1:2*j)=[2*idx(j)-1 2*idx(j)];
    end
end

%% --- Material properties ---
E  = 210e9;  % Young's modulus [Pa]
nu = 0.3;    % Poisson's ratio

%% --- Gauss points for 2x2 integration ---
K = zeros(ndof,ndof);
for i=1:num_ele
    Ke = elementalStiffnessPlate(coordinatesFEM(elementsFEM(i,:)',:),E,nu);
    idx = elementDofs(i,:);
    K(idx,idx)=K(idx,idx)+Ke;
end

% Apply Boundary Condition
K(dofs_bottom,:)=[];
K(:,dofs_bottom)=[];
F(dofs_bottom)=[];
mu_f(dofs_bottom)=[];

%% Setup of inverse problem
% forward operator
G = K\C';
G = G';

% "exact" posterior covariance
gamma_pos = gamma_prior_f-gamma_prior_f*G'*((G*gamma_prior_f*G'+gamma_obs)\G)*gamma_prior_f';

%% Approximations 
% LIS basis and singular values
[Omega,Delta,Phi]=svd((S_obs\G)*S_pr);
V = S_pr*Phi(:,1:m);
W = G'*((S_obs\Omega(:,1:m))/Delta(1:m,1:m));
delta = diag(Delta);

% actual errors
cov_error_ROM = zeros(1,m-1);
cov_error_OLR = zeros(1,m-1);
HessErr_ROM = zeros(1,m-1);
HessErr_OLR = zeros(1,m-1);
% error bound and constants
C1_ROM = zeros(1,m-1);
C1_OLR = zeros(1,m-1);
b_ROM = zeros(1,m-1);
t_OLR = zeros(1,m-1);
for r=1:m-1
    %% ROM approximation
    K_hat = W(:,1:r)'*K*V(:,1:r);
    G_hat = C*V(:,1:r)*(K_hat\W(:,1:r)');
    G_ROM_cell{r}=G_hat;
    gamma_pos_ROM{r} = gamma_prior_f-gamma_prior_f*G_hat'*((G_hat*gamma_prior_f*G_hat'+gamma_obs)\G_hat)*gamma_prior_f';
    cov_error_ROM(r)=norm(gamma_pos_ROM{r}-gamma_pos,2);
    
    %% OLR approximation 
    G_r = G*V(:,1:r)*W(:,1:r)';
    G_OLR_cell{r}=G_r;
    gamma_pos_OLR{r} = gamma_prior_f-gamma_prior_f*G_r'*((G_r*gamma_prior_f*G_r'+gamma_obs)\G_r)*gamma_prior_f';
    cov_error_OLR(r)=norm(gamma_pos-gamma_pos_OLR{r},2);
    
    %% Error constants
    G_hat=G_ROM_cell{r};
    HessErr_ROM(r)=norm((S_obs\(G-G_hat))*S_pr,2);
    
    preHess = (S_obs\G)*S_pr;
    preHessROM = (S_obs\G_hat)*S_pr;
    preHessOLR = (S_obs\G_OLR_cell{r})*S_pr;

    i=r;

    HessErr_OLR(r)= norm(preHess-preHessOLR,2);
    norm_hess = norm(preHess,2);
    norm_hess_ROM = norm(preHessROM,2);
    norm_hess_OLR = norm(preHessOLR,2);

    R = sqrt(Delta(1:i,1:i)'*Delta(1:i,1:i))*W(:,1:i)'*K*V(:,1:i);

    C1_ROM(i) = norm(S_pr,2)^2*(norm((eye(m)+preHess*preHess')\preHess,2)+ ...
    norm_hess_ROM*norm_hess*(norm_hess_ROM+norm_hess) ...
    +norm((eye(m)+preHessROM*preHessROM')\preHessROM,2));
    
    C1_OLR(i) = norm(S_pr,2)^2*(norm((eye(m)+preHess*preHess')\preHess,2)+ ...
    norm_hess_OLR*norm_hess*(norm_hess_OLR+norm_hess) ...
    +norm((eye(m)+preHessOLR*preHessOLR')\preHessOLR,2));
    
    t_OLR(i)=norm(sqrt(Delta(i+1:end,i+1:end)'*Delta(i+1:end,i+1:end)),2);
    b_ROM(i)=t_OLR(i)+...
            norm(Omega(:,i+1:end)*Omega(:,i+1:end)'*(S_obs\C)*(V(:,1:i)/R)*sqrt(Delta(1:i,1:i)'*Delta(1:i,1:i)),2)*norm(1,2);
       
end




%% Posterior mean
error_ROM=zeros(N_rep,m-1);
error_OLR=zeros(N_rep,m-1);
C2_ROM = zeros(1,m-1);
C2_OLR = zeros(1,m-1);

L_dagger = (S_pr'*S_pr)\S_pr';

for j=1:N_rep
    y = G*F+S_obs*randn(m,1);
    mu_full = mu_f+ gamma_pos*G'*(gamma_obs\(y-G*mu_f));
    for r=1:m-1
        G_hat = G_ROM_cell{r};
        mu_LI = mu_f+gamma_prior_f*G_hat'*((G*gamma_prior_f*G'+gamma_obs)\(y-G_hat*mu_f));
        G_hat = G_OLR_cell{r};
        mu_OLR = mu_f+gamma_prior_f*G_hat'*((G*gamma_prior_f*G'+gamma_obs)\(y-G_hat*mu_f));

        error_ROM(j,r) = norm(mu_full-mu_LI,2);
        
        error_OLR(j,r) = norm(mu_full-mu_OLR,2);

        %% Bounds
        i=r;
        C2_ROM(i)=C1_ROM(i)*norm(G'*(S_obs\(y-G*mu_f)),2)+...
            norm(gamma_pos_ROM{i}*G'/S_obs,2)*norm(L_dagger*mu_f,2)+...
            norm(gamma_pos_ROM{i}*L_dagger',2)*norm(S_obs\(y-G_ROM_cell{i}*mu_f),2);

        C2_OLR(i)=C1_OLR(i)*norm(G'*(S_obs\(y-G*mu_f)),2)+...
            norm(gamma_pos_OLR{i}*G'/S_obs,2)*norm(L_dagger*mu_f,2)+...
            norm(gamma_pos_OLR{i}*L_dagger',2)*norm(S_obs\(y-G_OLR_cell{i}*mu_f),2);
    end
end

mean_error_ROM = mean(error_ROM,1);
mean_error_OLR = mean(error_OLR,1);


%% PLOTS

width = 7.25;
height = 6;
alpha = 0.25;
ROM_color = (1-alpha)*[0.4660 0.6740 0.1880]+alpha*[1 1 1];
OLR_color = (1-alpha)*[0.8500 0.3250 0.0980]+alpha*[1 1 1];

figure;

semilogy(delta,"LineWidth",2)
set(gca,"FontSize",20)
box off
title("Singular values $\delta_i$","Interpreter","latex","FontSize",28)
xlabel("$i$","Interpreter","latex")
axis([1 6 1e-1 5*1e1])
%yticks([10^(-15) 10^(-10) 10^(-5) 10^0])

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0.5 0.5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);

%exportgraphics(gcf, 'deltaPlate.pdf', 'ContentType', 'vector');

fig = figure;
semilogy(b_ROM,"x","Color",ROM_color,"LineWidth",2,"MarkerSize",8)
set(gca,"FontSize",20)
box off
hold on
semilogy(t_OLR,"o","Color",OLR_color,"LineWidth",2)
semilogy(HessErr_ROM,"Color",ROM_color,"LineWidth",2)
semilogy(HessErr_OLR,"Color",OLR_color,"LineWidth",2)
axis([1 5.1 1e-1 5*1e1])
%yticks([10^(-15) 10^(-10) 10^(-5) 10^0])
title("$S_{\mathrm{obs}}^{-1}(G-\hat{G}(r))S_{\mathrm{pr}}$ error and bound","Interpreter","latex","FontSize",28)
legend("ROM bound","OLR bound","ROM","OLR","Interpreter","Latex","Location","northeast")
legend boxoff
xlabel("Approximation rank $r$","Interpreter","latex")

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0.5 0.5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);

%exportgraphics(gcf, 'hessPlate.pdf', 'ContentType', 'vector');

figure;
semilogy(C2_ROM.*b_ROM,"x","Color",ROM_color,"LineWidth",2,"MarkerSize",8)
set(gca,'FontSize',20)
box off
hold on
semilogy(C2_OLR.*t_OLR,"o","Color",OLR_color,"LineWidth",2)
semilogy(mean_error_ROM,'Color',ROM_color,'LineWidth',2)
semilogy(mean_error_OLR,'Color',OLR_color,'LineWidth',2)
legend('ROM bound','OLR bound','ROM','OLR','Location','northwest')
legend boxoff
title('$\mu_{\mathrm{pos}}$ error and bound','Interpreter','latex','FontSize',28)
axis([1 5 1e4 1e20])
xlabel('Approximation rank $r$','Interpreter','latex')

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0.5 0.5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);

%exportgraphics(gcf, 'meanPlate.pdf', 'ContentType', 'vector');

figure;
semilogy(C1_ROM.*b_ROM,"x","Color",ROM_color,"LineWidth",2,"MarkerSize",8)
set(gca,'FontSize',20)
box off
hold on
semilogy(C1_OLR.*t_OLR,"o","Color",OLR_color,"LineWidth",2)
semilogy(cov_error_ROM,'Color',ROM_color,'LineWidth',2)
semilogy(cov_error_OLR,'Color',OLR_color,'LineWidth',2)
title('$\Gamma_{\mathrm{pos}}$ error and bound for $p=\infty$','Interpreter','latex','FontSize',28)
xlabel('Approximation rank $r$','Interpreter','latex')
legend('ROM bound','OLR bound','ROM','OLR','Location','southwest')
legend boxoff
axis([1 5 1e4 1e20])

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0.5 0.5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);

%exportgraphics(gcf, 'covPlate.pdf', 'ContentType', 'vector');
