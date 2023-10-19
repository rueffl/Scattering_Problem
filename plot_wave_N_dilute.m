%% Set structure setting
clear all
format long

% Settings for the material's structure
k_tr = 0; % truncation parameters as in remark 3.3
N = 2; % number of the resonator
spacing = 1000; lij = ones(1,N-1).*spacing; % spacing between the resonators
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
uin = @(x,t,n) exp((k0).*x+w0.*t).*(x<xm(1)).*(n==0); % incident wave
vin = @(x) exp(k0*x);
G = @(k,x) exp(sqrt(-1)*k*abs(x))./(2*sqrt(-1)*k); % Green's function

% Define evaluation points
len_xs = 800;
len_zs = 20;
xs = linspace(xm(1)-1,xp(end)+1,len_xs);
zs = zeros(N,len_zs);
for i = 1:N
    zs(i,:) = linspace(xm(i),xp(i),len_zs);
end

%% Calculate the scattered wave field taking the left and right incident wave fields into account 

uin_r = @(x,t,n) 0; vin_r = @(x,n) 0; dx_vin_r = @(x,n) 0; % zero right incident field

alpha = zeros(N,2*k_tr+1); beta = zeros(N,2*k_tr+1);

for i = 1:N

    zi = z(i);
    C = getC(k_tr, i, w_op, Omega, rs, ks, vr);
    [f_ni,lambdas] = eig(C,'vector');
    lambdas = flip(sqrt(lambdas));

    % Define relevant functions for the incident wave field
    if i == 1
        uin_l = @(x,t,n) exp(sqrt(-1)*((k0).*x+w0.*t)).*(x<xm(1)).*(n==0); % left incident wave
        vin_l = @(x,n) exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % n-th mode of left incident wave
        dx_vin_l = @(x,n) sqrt(-1)*k0*exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % derivative of left incident wave
    else
        alpha_i_old = alpha_i;
        uin_l = @(x,t,n) alpha_i_old(n+k_tr+1)*exp(sqrt(-1)*x*(w_op+n*Omega)/v0).*exp(sqrt(-1)*n*Omega*t).*(xp(i-1)<x && x<xm(i)); % left incident wave comming from the previous resonator
        vin_l = @(x,n) alpha_i_old(n+k_tr+1)*exp(sqrt(-1)*x*(w_op+n*Omega)/v0).*(xp(i-1)<x && x<xm(i)); % n-th mode of the left incident wave
        dx_vin_l = @(x,n) alpha_i_old(n+k_tr+1)*sqrt(-1)*(w_op+n*Omega)/v0*exp(sqrt(-1)*x*(w_op+n*Omega)/v0).*(xp(i-1)<x && x<xm(i)); % derivative of left incident wave
        uin = @(x,t,n) uin(x,t,n) + uin_l(x,t,n);
    end
    
    % Compute solution coefficients
    MatcalA = getMatcalA(1,[],xm(i),xp(i),k_tr,w_op,Omega,rs(i,:),ks(i,:),vr,delta,v0); % matrix \mathcal{A}
    MatcalF = getF_lr(k_tr, 1, delta, xm(i), xp(i), dx_vin_l, dx_vin_r); % vector \mathcal{F}
    sol = MatcalA\MatcalF; % solve for the interior coefficients, vector \mathbf{w}
    as = flip(sol(1:2:end)); bs = flip(sol(2:2:end)); 

    % Calculate alpha_n^i and beta_n^i % for every n = -K, ..., K
    alpha_i = zeros(1,2*k_tr+1); beta_i = zeros(1,2*k_tr+1); % prepare for values of alpha and beta for a fixed resonator D_i for all n
    for n = -k_tr:k_tr % iterate over the modes

        sum_ni = 0;
        kn = (w_op+n*Omega)/v0; % n-th wave number outside of the resonators
        for j = -k_tr:k_tr
            sum_ni = sum_ni + (as(j+k_tr+1)*exp(sqrt(-1)*lambdas(j+k_tr+1)*zi)+bs(j+k_tr+1)*exp(-sqrt(-1)*lambdas(j+k_tr+1)*zi))*f_ni(j+k_tr+1);
        end
        alpha_ni = sum_ni*exp(-sqrt(-1)*kn*zi);
        beta_ni = sum_ni - vin_l(zi-0.0000001,n);
        beta_ni = beta_ni*exp(sqrt(-1)*kn*zi);
    
        alpha_i(n+k_tr+1) = alpha_ni;
        beta_i(n+k_tr+1) = beta_ni;
    
    end

    alpha(i,:) = alpha_i;
    beta(i,:) = beta_i;

end

% Define the scattered wave field
u_sc = @(x,t) get_usc(x,t,alpha,beta,vin,z,w_op,Omega,v0,k_tr,uin);

for i = 1:(N+1)

    % Define evaluation points and evaluate the scattered wave at them
    len_xs = 800;
    if i == 1
        xs = linspace(z(i)-4,z(i),len_xs);
    elseif i == N+1
        xs = linspace(z(i-1),z(i-1)+4,len_xs);
    else
        xs = linspace(z(i-1),z(i),len_xs);
    end
    u_x = zeros(1,len_xs);
    for k = 1:len_xs
        u_x(k) = u_sc(xs(k),0);
    end

    % Plot the evaluated scattered wave field
    plot(xs,u_x,'.')
    hold on

end



%% Functions used in this file

function [usc_x] = get_usc(x,t,alphas,betas,vin,z,w_op,Omega,v0,k_tr,uin)
%GET_VN  Evaluates v_n at x for a given n
%   x:      spatial coordinate
%   t:      temporal coordinate
%   alphas: coefficients \alpha_n^i, for all i=1,...,N and n=-K,...,K
%   betas:  coefficients \beta_n^i, for all i=1,...,N and n=-K,...,K
%   vin:    incidnt wave field to the system
%   z:      location of resonators
%   w_op:   operating frequency 
%   Omega:  frequency of the time-modulated material parameters
%   v0:     wave speed outside of the resonators
%   k_tr:   truncation parameter

    usc_x = 0;
    for n = -k_tr:k_tr
        vn_sc_x = get_vn(x,alphas(:,n+k_tr+1),betas(:,n+k_tr+1),(w_op+n*Omega)/v0,vin,z);
        usc_x = usc_x + vn_sc_x*exp(sqrt(-1)*(w_op+n*Omega)*t)+uin(x,t,n);
    end

end

function [vn_x] = get_vn(x,alphas_n,betas_n,kn,vin,z)
%GET_VN  Evaluates v_n at x for a given n
%   x:          spatial coordinate
%   alphas_n:   coefficients \alpha_n^i, for all i=1,...,N
%   betas_n:    coefficients \beta_n^i, for all i=1,...,N
%   kn:         n-th wave number outside of D
%   vin:        incidnt wave field to the system
%   z:          location of resonators

    if x<=z(1)
        vn_x = betas_n(1)*exp(-sqrt(-1)*kn*x)+vin(x);
    elseif x>z(end)
        vn_x = alphas_n(end)*exp(sqrt(-1)*kn*x);
    else
        for i = 1:(length(betas_n)-1)
            if z(i)<x && x<=z(i+1)
                vn_x = alphas_n(i)*exp(sqrt(-1)*kn*x)+betas_n(i+1)*exp(-sqrt(-1)*kn*x);
            end
        end
    end

end

