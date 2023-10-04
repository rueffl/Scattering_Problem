%% Set structure setting
clear all
format long

% Settings for the material's structure
k_tr = 4; % truncation parameters as in remark 3.3
k_tr_n = k_tr;
k_tr_m = 0;
N = 3; % number of the resonator
spacing = 10; lij = ones(1,N-1).*spacing; % spacing between the resonators
len = 0.0005; li = ones(1,N).*len; % length of the resonator
L = sum(li)+sum(lij); % length of the unit cell
Ls = zeros(2*N-1,1);
Ls(1:2:end) = li;
Ls(2:2:end) = lij;
xipm = [0,cumsum(Ls)']+1; % all boundary points
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
% C = make_capacitance_finite(N,lij); % capacitance matrix
% w_muller = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0,lij,xm,xp); % subwavelength resonant frequencies
% w_res = w_muller(real(w_muller)>=0); % positive subwavelength resonant frequencies
w_res = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0); % non-zero subwavelength resonant frequency
w_op = w_res(1)+ 0.0002; % operating frequency
w0 = w_res(1); % quasifrequency of incident wave
% w0 = w_op;
k_op = w_op/v0; % operating wave number outside of the resonator
k0 = w0/v0; % wave number of incident frequency

% Define relevant functions
uin = @(x,t) exp((k0).*x+w0.*t); % incident wave
G = @(k,x) exp(sqrt(-1)*k*abs(x))./(2*sqrt(-1)*k); % Green's function

% Define evaluation points
len_xs = 800;
len_zs = 80;
xs = linspace(xm(1)-1,xp(end)+1,len_xs);
zs = zeros(N,len_zs);
for i = 1:N
    zs(i,:) = linspace(xm(i),xp(i),len_zs);
end

%% Calculate the scattered wave field 

Lambdan = zeros(N,2*k_tr+1);
fig = figure();
hold on
xs = linspace(xm(1)-spacing,xp(1),len_xs);
plot(xs,uin(xs,t),'-')

for i = 1:N
    
    if i == 1
        MatcalF = getF(k_tr, 1, delta, k_op, k0, xm(1)); % vector \mathcal{F}
    else
        MatcalF = getFi(k_tr, delta, w_res, w_op, Lambdan_xneg, z(i-1), xm(i), Omega, v0); % vector \mathcal{F}_i
    end
    MatcalA = getMatcalA(1, [], xm(i), xp(i), k_tr, w_op, Omega, rs(i,:), ks(i,:), vr, delta, v0); % matrix \mathcal{A}
    sol = MatcalA\MatcalF; % solve for the interior coefficients, vector \mathbf{w}
    if i > 1
        Lambdan_old = @(x) Lambdan(x);
    end
    Lambdan = @(x) get_Lambdas_N1(z,xm(i),xp(i),k_tr,w_op,Omega,rs,ks,vr,v0,delta,k_op,k0,x,i,sol); % frequency-scattering coefficients
    Lambdan_xneg = get_Lambdas_N1(z,xm(i),xp(i),k_tr,w_op,Omega,rs,ks,vr,v0,delta,k_op,k0,xm(i)-0.01,i,sol);

    if i == 1
        u_tot = @(x,t) get_usc_i(x, t, Lambdan, k_tr, w_op, w_res, Omega, v0, z(1));% + uin(x,t);
    else
        u_tot_old = @(x,t) u_tot(x,t);%get_usc_i(x, t, Lambdan_old, k_tr, w_op, w_res, Omega, v0, z(i-1));
        u_tot = @(x,t) get_usc_i(x, t, Lambdan, k_tr, w_op, w_res, Omega, v0, z(i)) + (x~=z(i)).*u_tot_old(x,t);
    end
    
    if i < N
        xs = linspace(xp(i),xp(i+1),len_xs);
%         xs = linspace(xm(i),xm(i+1),len_xs);
    else
        xs = linspace(xp(i),xp(i)+spacing+len,len_xs);
%         xs = linspace(xm(i),xm(i+1)+spacing+len,len_xs);
    end
    ux = zeros(1,len_xs);
    for j = 1:len_xs
        ux(j) = u_tot(xs(j),t);
    end

    plot(xs,ux,'b-',LineWidth=2)

end





%% Functions

function usc = get_usc_i(x, t, Lambdan, k_tr, w_op, w_res, Omega, v0, zi)
% GET_USC_I  Calculates the scattered wave u^{sc}_i at x and time t
%   x:          spatial coordinate
%   t:          temporal coordinate
%   Lambdan:    frequency-scattering coefficient
%   k_tr:       truncation parameter
%   w_op:       operating wave frequency
%   w_res:      resonant frequency
%   Omega:      frequency of time-modulation
%   v0:         wave speed outside of the resonators
%   zi:         centre of resonator

    G = @(k,x) exp(sqrt(-1)*k*abs(x))./(2*sqrt(-1)*k); % Green's function
    usc = 0;
    select = @(M,idx) M(idx);
    for n = -k_tr:k_tr
        kn = (w_op+n*Omega)/v0;
        Lambda = @(x) select(Lambdan(x),n+k_tr+1);
        usc = usc + Lambda(x)*G(kn,x-zi)*exp(sqrt(-1)*(w_res+n*Omega)*t);
%         usc = usc + Lambdan(n+k_tr+1)*G(kn,x-zi)*exp(sqrt(-1)*(w_res+n*Omega)*t);
    end
end







