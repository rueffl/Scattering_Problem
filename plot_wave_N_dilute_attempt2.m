%% Set structure setting
clear all
format long

% Settings for the material's structure
k_tr = 3; % truncation parameters as in remark 3.3
N = 2; % number of the resonator
spacing = 100; lij = ones(1,N-1).*spacing; % spacing between the resonators
len = 0.05; li = ones(1,N).*len; % length of the resonator
L = sum(li)+sum(lij); % length of the unit cell
Ls = zeros(2*N-1,1);
Ls(1:2:end) = li;
Ls(2:2:end) = lij;
xipm = [0,cumsum(Ls)']-(len*(N/2)+spacing*(N-1)/2); % all boundary points, make sure the resonators are aligned symmetrically wrt 0
xm = xipm(1:2:end); % LHS boundary points
xp = xipm(2:2:end); % RHS boundary points
z = (xm+xp)./2; % centers of resonators
delta = 0.0001; % small contrast parameter
t = 0; % time

vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators

% Setting for the material's time-modulation
Omega = 0.03; % modulation frequency
T = 2*pi/Omega;
phase_kappa = zeros(1,N); % modulation phases of kappa
phase_rho = zeros(1,N); % modulation phases of rho
for i = 1:(N-1)
    phase_kappa(i+1) = pi/i;
    phase_rho(i+1) = pi/i;
end
epsilon_kappa = 0; % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho
rs = []; % Fourier coefficients of 1/rho
ks = []; % Fourier coefficients of 1/kappa
for j = 1:N
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
    ks = [ks; ks_j];
    rs = [rs; rs_j];
end

% Calculate subwavelength resonant frequency
w_res = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0); % non-zero subwavelength resonant frequency
w_op = w_res(1)+ 0.0002; % operating frequency
w0 = w_res(1); % quasifrequency of incident wave
% w0 = w_op;
k_op = w_op/v0; % operating wave number outside of the resonator
k0 = w0/v0; % wave number of incident frequency

% Define relevant functions
G = @(k,x) exp(sqrt(-1)*k*abs(x))./(2*sqrt(-1)*k); % Green's function



%% Add solution for left incident wave and right incident wave

all_alphas_l = zeros(N,2*k_tr+1); all_betas_l = zeros(N,2*k_tr+1);
all_alphas_r = zeros(N,2*k_tr+1); all_betas_r = zeros(N,2*k_tr+1);

% left incident wave
for i = 1:N

    % left incident wave
    if i == 1
        vin_l = @(x,n) exp(sqrt(-1)*(k0).*x).*(x<xm(i)).*(n==0); % n-th mode of left incident wave
        dx_vin_l = @(x,n) sqrt(-1)*k0*exp(sqrt(-1)*(k0).*x).*(x<xm(i)).*(n==0); % derivative of left incident wave
%         vin_l = @(x,n) 0; dx_vin_l = @(x,n) 0;
        uin_r = @(x,t,n) 0; vin_r = @(x,n) 0; dx_vin_r = @(x,n) 0;
    else
        uin_l = @(x,t,n) alphas(n+k_tr+1)*exp(sqrt(-1)*((w_op+n*Omega)/v0*x+(n*Omega+w_op)*t)); % left incident wave from D_{i-1}
        vin_l = @(x,n) alphas(n+k_tr+1)*exp(sqrt(-1)*((w_op+n*Omega)/v0*x)); % n-th mode of left incident wave
        dx_vin_l = @(x,n) sqrt(-1)*(w_op+n*Omega)/v0*alphas(n+k_tr+1)*exp(sqrt(-1)*((w_op+n*Omega)/v0*x)); % derivative of n-th mode of left incident wave
    end
    
    % Compute solution coefficients
    MatcalA = getMatcalA(1,[],xm(i),xp(i),k_tr,w_op,Omega,rs(i,:),ks(i,:),vr,delta,v0); % matrix \mathcal{A}
    MatcalF = getF_lr(k_tr, 1, delta, xm(i), xp(i), dx_vin_l, dx_vin_r); % vector \mathcal{F}
    sol_l = MatcalA\MatcalF; % solve for the interior coefficients, vector \mathbf{w}

    % Compute the eigenvectors and the sqrt of the eigenvalues of C_i
    C = getC(k_tr, i, w_op, Omega, rs, ks, vr);
    [fs,lambdas] = eig(C,'vector');
    lambdas = sqrt(lambdas);

    [alphas,betas] = get_alpha_beta(sol_l,lambdas,fs,xm(i),xp(i),w_op,Omega,v0,k_tr,vin_r,vin_l); % calculate exterior coefficients alpha and beta for D_i 
    all_alphas_l(i,:) = alphas; all_betas_l(i,:) = betas;

end

len_xs = 8000;
% all_xs = zeros(1,(N+1)*len_xs);
all_xs = linspace(xm(1)-4,xp(end)+4,len_xs);
all_uxs_l = zeros(1,len_xs); all_uxs_r = zeros(1,len_xs);
all_uxs_in_l = zeros(1,len_xs); all_uxs_in_r = zeros(1,len_xs);
t = 0;

% evaluate u_sc and u_in between each resonators
for i = 1:(N+1)

    if i == 1
%         all_xs(1:len_xs) = linspace(xm(1)-4,xm(1),len_xs); % domain left of D_1
        v_sc = @(x,n) all_betas_l(1,n+k_tr+1)*exp(-sqrt(-1)*x*(w_op+n*Omega)/v0);
        u_sc = @(x,t) get_u_from_v(x,t,v_sc,k_tr,w_op,Omega); % scattered wave field left of D_1
        u_in = @(x,t) exp(sqrt(-1)*((k0).*x+w0.*t)).*(x<=xm(1)); % left incident wave of D_1
    elseif i < N+1
%         all_xs((i-1)*len_xs+1:i*len_xs) = linspace(xp(i-1),xm(i),len_xs); % domain between D_{i-1} and D_i
        v_sc = @(x,n) all_alphas_l(i-1,n+k_tr+1)*exp(sqrt(-1)*x*(w_op+n*Omega)/v0)+all_betas_l(i,n+k_tr+1)*exp(-sqrt(-1)*x*(w_op+n*Omega)/v0);
        u_sc = @(x,t) get_u_from_v(x,t,v_sc,k_tr,w_op,Omega); % scattered wave field between D_{i-1} and D_i
        v_in = @(x,n) all_alphas_l(i-1,n+k_tr+1)*exp(sqrt(-1)*x*(w_op+n*Omega)/v0);
        u_in = @(x,t) get_u_from_v(x,t,v_in,k_tr,w_op,Omega).*(xp(i-1)<=x && x<=xm(i)); % left incident wave field of D_{i-1}
    else
%         all_xs(N*len_xs+1:(N+1)*len_xs) = linspace(xp(N),xp(N)+4,len_xs); % domain right of D_N
        v_sc = @(x,n) all_alphas_l(N,n+k_tr+1)*exp(sqrt(-1)*x*(w_op+n*Omega)/v0);
        u_sc = @(x,t) get_u_from_v(x,t,v_sc,k_tr,w_op,Omega); % scattered wave field right of D_N
        v_in = @(x,n) 0;%all_alphas_l(N,n+k_tr+1)*exp(sqrt(-1)*x*(w_op+n*Omega)/v0);
        u_in = @(x,t) 0;%get_u_from_v(x,t,v_in,k_tr,w_op,Omega); % left incident wave field of D_N
    end

    for k = 1:len_xs
        x = all_xs(k);
        all_uxs_l(k) = all_uxs_l(k) + u_sc(x,t);
        all_uxs_in_l(k) = all_uxs_in_l(k) + u_in(x,t);
    end
    
end
% 
% 
% % right incident wave
% for i = N:(-1):1
% 
%     % right incident wave
%     if i == N
%         uin_l = @(x,t,n) 0; vin_l = @(x,n) 0; dx_vin_l = @(x,n) 0;
%         uin_r = @(x,t,n) 0; vin_r = @(x,n) 0; dx_vin_r = @(x,n) 0; % no right incident wave to the system
%         uin_r = @(x,t,n) exp(sqrt(-1)*(-(k0).*x+w0.*t)).*(x>xp(end)).*(n==0); % right incident wave
%         vin_r = @(x,n) exp(sqrt(-1)*(-(k0).*x)).*(x>xp(end)).*(n==0); % n-th mode of right incident wave
%         dx_vin_r = @(x,n) -sqrt(-1)*k0*exp(sqrt(-1)*(-(k0).*x)).*(x>xp(end)).*(n==0); % derivative of right incident wave
%     else
%         uin_r = @(x,t,n) betas(n+k_tr+1)*exp(-sqrt(-1)*((w_op+n*Omega)/v0*x+(n*Omega+w_op)*t)); % right incident wave from D_{i-1}
%         vin_r = @(x,n) betas(n+k_tr+1)*exp(-sqrt(-1)*((w_op+n*Omega)/v0*x)); % n-th mode of right incident wave
%         dx_vin_r = @(x,n) -sqrt(-1)*(w_op+n*Omega)/v0*betas(n+k_tr+1)*exp(-sqrt(-1)*((w_op+n*Omega)/v0*x)); % derivative of n-th mode of right incident wave
%     end
%     
%     % Compute solution coefficients
%     MatcalA = getMatcalA(1,[],xm(i),xp(i),k_tr,w_op,Omega,rs(i,:),ks(i,:),vr,delta,v0); % matrix \mathcal{A}
%     MatcalF = getF_lr(k_tr, 1, delta, xm(i), xp(i), dx_vin_l, dx_vin_r); % vector \mathcal{F}
%     sol_r = MatcalA\MatcalF; % solve for the interior coefficients, vector \mathbf{w}
% 
%     % Compute the eigenvectors and the sqrt of the eigenvalues of C_i
%     C = getC(k_tr, i, w_op, Omega, rs, ks, vr);
%     [fs,lambdas] = eig(C,'vector');
%     lambdas = sqrt(lambdas);
% 
%     [alphas,betas] = get_alpha_beta(sol_r,lambdas,fs,xm(i),xp(i),w_op,Omega,v0,k_tr,vin_r,vin_l); % calculate exterior coefficients alpha and beta for D_i 
%     all_alphas_r(i,:) = alphas; all_betas_r(i,:) = betas;
% 
% end

% % evaluate u_sc and u_in between each resonators
% for i = 1:(N+1)
% 
%     if i == 1 % domain left of D_1
%         v_sc = @(x,n) all_betas_r(1,n+k_tr+1)*exp(-sqrt(-1)*x*(w_op+n*Omega)/v0);
%         u_sc = @(x,t) get_u_from_v(x,t,v_sc,k_tr,w_op,Omega); % scattered wave field left of D_1
%         v_in = @(x,n) all_betas_r(1,n+k_tr+1)*exp(-sqrt(-1)*x*(w_op+n*Omega)/v0); % n-th mode of right incident wave of D_1
%         u_in = @(x,t) get_u_from_v(x,t,v_in,k_tr,w_op,Omega); % right incident wave of D_1
%     elseif i < N+1 % domain between D_{i-1} and D_i
%         v_sc = @(x,n) all_alphas_r(i-1,n+k_tr+1)*exp(sqrt(-1)*x*(w_op+n*Omega)/v0)+all_betas_r(i,n+k_tr+1)*exp(-sqrt(-1)*x*(w_op+n*Omega)/v0);
%         u_sc = @(x,t) get_u_from_v(x,t,v_sc,k_tr,w_op,Omega); % scattered wave field between D_{i-1} and D_i
%         v_in = @(x,n) all_betas_r(i,n+k_tr+1)*exp(-sqrt(-1)*x*(w_op+n*Omega)/v0);
%         u_in = @(x,t) get_u_from_v(x,t,v_in,k_tr,w_op,Omega); % left incident wave field of D_{i-1}
%     else % domain right of D_N
%         v_sc = @(x,n) all_alphas_r(N,n+k_tr+1)*exp(sqrt(-1)*x*(w_op+n*Omega)/v0);
%         u_sc = @(x,t) get_u_from_v(x,t,v_sc,k_tr,w_op,Omega); % scattered wave field right of D_N
%         u_in = @(x,t) 0; % right incident wave field of D_N
%     end
% 
%     for k = 1:len_xs
%         x = all_xs((i-1)*len_xs+k);
%         all_uxs_r((i-1)*len_xs+k) = u_sc(x,t);
%         all_uxs_in_r((i-1)*len_xs+k) = u_in(x,t);
%     end
%     
% end

all_uxs_sc = all_uxs_l + all_uxs_r;
all_uxs_in = all_uxs_in_l + all_uxs_in_r;
all_uxs = zeros(1,length(all_uxs_in));
for k = 1:length(all_uxs_in)
    all_uxs(k) = all_uxs_sc(k)+all_uxs_in(k);
end

% Create plots of the scattered, incident and total wave field (left)
[fig1, fig2, fig3] = create_plots(all_xs,t,all_uxs_l,all_uxs_in_l,all_uxs_l+all_uxs_in_l);

% % Create plots of the scattered, incident and total wave field (right)
% [fig4, fig5, fig6] = create_plots(all_xs,t,all_uxs_r,all_uxs_in_r,all_uxs_r+all_uxs_in_r);
% 
% % Create plots of the scattered, incident and total wave field
% [fig7, fig8, fig9] = create_plots(all_xs,t,all_uxs_sc,all_uxs_in,all_uxs);

%% functions

function [u_xt] = get_u_from_v(x,t,vn,k_tr,w_op,Omega)
%GET_U_FROM_V  Adds up all of the modes to obtain u
%   x:      spatial coordinate
%   t:      temporal coordinate
%   vn:     modes of incident wave field
%   k_tr:   truncation parameter
%   w_op:   operating frequency 
%   Omega:  frequency of the time-modulated material parameters 

    u_xt = 0;
    for n = -k_tr:k_tr
        u_xt_n = vn(x,n)*exp(sqrt(-1)*(w_op+n*Omega)*t);
        u_xt = u_xt + u_xt_n;
    end

end

function [fig1, fig2, fig3] = create_plots(all_xs,t,all_uxs_sc,all_uxs_in,all_uxs)
%CREATE_PLOTS  Creates the plots of the incident, scattered and total wave fields
%   all_xs:         spatial evaluation points
%   t:              temporal coordinate
%   all_uxs_sc:     scattered wave field
%   all_uxs_in:     incident wave field
%   all_uxs:        total wave field
    
    % create plot of scattered wave field
    fig1 = figure();
    fig1.Position = [263,725,982,352];
    subplot(1,2,1)
    set(gca,'FontSize',14)
    hold on
    plot(all_xs,real(all_uxs_sc),'-','DisplayName','Exact',markersize=8,linewidth=2)
    xlim([all_xs(1) all_xs(end)])
    xlabel('$x$',Interpreter='latex',FontSize=18)
    ylabel(strcat('Re$(u^{\mathrm{sc}}(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)
    
    subplot(1,2,2)
    set(gca,'FontSize',14)
    hold on
    plot(all_xs,imag(all_uxs_sc),'-','DisplayName','Exact',markersize=8,linewidth=2)
    xlim([all_xs(1) all_xs(end)])
    xlabel('$x$',Interpreter='latex',FontSize=18)
    ylabel(strcat('Im$(u^{\mathrm{sc}}(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)
    
    % create plot of incident wave field
    fig2 = figure();
    fig2.Position = [263,725,982,352];
    subplot(1,2,1)
    set(gca,'FontSize',14)
    hold on
    plot(all_xs,real(all_uxs_in),'-','DisplayName','Exact',markersize=8,linewidth=2)
    xlim([all_xs(1) all_xs(end)])
    xlabel('$x$',Interpreter='latex',FontSize=18)
    ylabel(strcat('Re$(u^{\mathrm{in}}(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)
    
    subplot(1,2,2)
    set(gca,'FontSize',14)
    hold on
    plot(all_xs,imag(all_uxs_in),'-','DisplayName','Exact',markersize=8,linewidth=2)
    xlim([all_xs(1) all_xs(end)])
    xlabel('$x$',Interpreter='latex',FontSize=18)
    ylabel(strcat('Im$(u^{\mathrm{in}}(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)
    
    % create plot of total wave field
    fig3 = figure();
    fig3.Position = [263,725,982,352];
    subplot(1,2,1)
    set(gca,'FontSize',14)
    hold on
    plot(all_xs,real(all_uxs),'-','DisplayName','Exact',markersize=8,linewidth=2)
    xlim([all_xs(1) all_xs(end)])
    xlabel('$x$',Interpreter='latex',FontSize=18)
    ylabel(strcat('Re$(u(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)
    
    subplot(1,2,2)
    set(gca,'FontSize',14)
    hold on
    plot(all_xs,imag(all_uxs),'-','DisplayName','Exact',markersize=8,linewidth=2)
    xlim([all_xs(1) all_xs(end)])
    xlabel('$x$',Interpreter='latex',FontSize=18)
    ylabel(strcat('Im$(u(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)

end




