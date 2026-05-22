clearvars

addpath src\
%% Parameters of the problem
rng(25);

%% Setup of the problem
Parameters({'beam_bool'},{false});
load('Parameters.mat')

x_dofs=applyBoundaryCondition(x_dofs,BC_dofs,'Coordinate');

LIP_Setup({'gamma_obs'},{0.001^2});
load('LIP_Setup.mat')

%% Generalized eigenvectors / LIS basis
S_obs=chol(gamma_obs,"lower");
[Omega,Delta,Phi]=svd((S_obs\G)*S_pr);
delta = diag(Delta);
[~,V,W]=calculateLISBasis();

%% Mean samples
N_rep = 1;

% Generate sample
[~,state_sample]=gen_samples(N_rep);

d_f_ROM=zeros(1,m-1);
d_f_OLR=zeros(1,m-1);

p_ROM =zeros(1,m-1);
p_OLR = zeros(1,m-1);

C1_ROM = zeros(1,m-1);
C1_OLR = zeros(1,m-1);
HessErr_ROM = zeros(1,m-1);
HessErr_OLR = zeros(1,m-1);
b_ROM = zeros(1,m-1);
t_OLR = zeros(1,m-1);

for i=1:m-1

    %% ROM reduced operator
    [gamma_pos_ROM{i},G_ROM_cell{i},d_f_ROM(i)]=solveReducedModel(V(:,1:i),W(:,1:i)); 
    p_ROM(i)=norm(gamma_pos_ROM{i}-gamma_pos,2);
    

    %% OLR Reduction
    [gamma_pos_OLR{i},G_OLR_cell{i},d_f_OLR(i)]=solveOLRA(V(:,1:i),W(:,1:i));
    p_OLR(i)=norm(gamma_pos_OLR{i}-gamma_pos,2);

    %% Error constants
    G_hat=G_ROM_cell{i};
    HessErr_ROM(i)=norm((S_obs\(G-G_hat))*S_pr,2);
    
    preHess = (S_obs\G)*S_pr;
    preHessROM = (S_obs\G_hat)*S_pr;
    preHessOLR = (S_obs\G_OLR_cell{i})*S_pr;

    HessErr_OLR(i)= norm(preHess-preHessOLR,2);
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
        norm(Omega(:,i+1:end)*Omega(:,i+1:end)'*(S_obs\C)*(V(:,1:i)/R)*sqrt(Delta(1:i,1:i)'*Delta(1:i,1:i))*Phi(:,1:i)',2)*norm(1,2);
end

%% Posterior Mean Analysis

error_ROM = zeros(N_rep,m-1);
error_OLR = zeros(N_rep,m-1);

C2_ROM = zeros(1,m-1);
C2_OLR = zeros(1,m-1);

L_dagger = S_pr'/(S_pr*S_pr');

for j=1:N_rep
    ysam = C*state_sample(:,j)+sqrt(gamma_obs)*randn(m,1);
    % Full model mean
    mu_full = meanCalculation(G,ysam);

    for i = 1:m-1
        % ROM approximation
        mu_ROM = meanCalculation(G_ROM_cell{i},ysam);


        % OLR approximation
        mu_OLR = meanCalculation(G_OLR_cell{i},ysam);

        error_ROM(j,i) = norm(mu_full-mu_ROM);
        error_OLR(j,i) = norm(mu_full-mu_OLR);

        %% Bounds
        C2_ROM(i)=C1_ROM(i)*norm(G'*(S_obs\(ysam-G*mu_f)),2)+...
            norm(gamma_pos_ROM{i}*G'/S_obs,2)*norm(L_dagger*mu_f,2)+...
            norm(gamma_pos_ROM{i}*L_dagger',2)*norm(S_obs\(ysam-G_ROM_cell{i}*mu_f),2);

        C2_OLR(i)=C1_OLR(i)*norm(G'*(S_obs\(ysam-G*mu_f)),2)+...
            norm(gamma_pos_OLR{i}*G'/S_obs,2)*norm(L_dagger*mu_f,2)+...
            norm(gamma_pos_OLR{i}*L_dagger',2)*norm(S_obs\(ysam-G_OLR_cell{i}*mu_f),2);

    end
    
end

mean_ROM = mean(error_ROM,1);

mean_OLR = mean(error_OLR,1);


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
axis([1 10 1e-4 1e1])
yticks([10^(-4) 10^(-3) 10^(-2) 10^(-1) 10^0])

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0.5 0.5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);

%exportgraphics(gcf, 'deltaBar.pdf', 'ContentType', 'vector');

figure;
semilogy(b_ROM,"x","Color",ROM_color,"LineWidth",2,"MarkerSize",8)
set(gca,"FontSize",20)
box off
hold on
semilogy(t_OLR,"o","Color",OLR_color,"LineWidth",2)
semilogy(HessErr_ROM,"Color",ROM_color,"LineWidth",2)
semilogy(HessErr_OLR,"Color",OLR_color,"LineWidth",2)
axis([1 9 1e-4 1e1])
yticks([10^(-4) 10^(-3) 10^(-2) 10^(-1) 10^0])
title("$S_{\mathrm{obs}}^{-1}(G-\hat{G}(r))S_{\mathrm{pr}}$ error and bound","Interpreter","latex","FontSize",28)
legend("ROM bound","OLR bound","ROM","OLR","Interpreter","Latex","Location","southwest")
legend boxoff
xlabel("Approximation rank $r$","Interpreter","latex")


set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0.5 0.5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);

%exportgraphics(gcf, 'hessBar.pdf', 'ContentType', 'vector');

figure;
semilogy(C2_ROM.*b_ROM,"x","Color",ROM_color,"LineWidth",2,"MarkerSize",8)
set(gca,'FontSize',20)
box off
hold on
semilogy(C2_OLR.*t_OLR,"o","Color",OLR_color,"LineWidth",2)
semilogy(mean_ROM,'Color',ROM_color,'LineWidth',2)
semilogy(mean_OLR,'Color',OLR_color,'LineWidth',2)
legend('ROM bound','OLR bound','ROM','OLR','Location','northwest')
legend boxoff
title('$\mu_{\mathrm{pos}}$ error and bound','Interpreter','latex','FontSize',28)
axis([1 9 1e0 1e15])
xlabel('Approximation rank $r$','Interpreter','latex')

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0.5 0.5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);

%exportgraphics(gcf, 'meanBar.pdf', 'ContentType', 'vector');

figure;
semilogy(C1_ROM.*b_ROM,"x","Color",ROM_color,"LineWidth",2,"MarkerSize",8)
set(gca,'FontSize',20)
box off
hold on
semilogy(C1_OLR.*t_OLR,"o","Color",OLR_color,"LineWidth",2)
semilogy(p_ROM,'Color',ROM_color,'LineWidth',2)
semilogy(p_OLR,'Color',OLR_color,'LineWidth',2)
title('$\Gamma_{\mathrm{pos}}$ error and bound for $p=\infty$','Interpreter','latex','FontSize',28)
xlabel('Approximation rank $r$','Interpreter','latex')
legend('ROM bound','OLR bound','ROM','OLR','Location','southwest')
legend boxoff
axis([1 9 1e0 1e15])

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0.5 0.5 width height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width height]);

%exportgraphics(gcf, 'covBar.pdf', 'ContentType', 'vector');
