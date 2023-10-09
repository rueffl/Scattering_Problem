%% Set structure setting
clear all
format long

% Settings for the material's structure
k_tr = 4; % truncation parameters as in remark 3.3
k_tr_n = k_tr;
k_tr_m = 0;
N = 6; % number of the resonator
spacing = 10; lij = ones(1,N-1).*spacing; % spacing between the resonators
len = 2; li = ones(1,N).*len; % length of the resonator
L = sum(li)+sum(lij); % length of the unit cell
Ls = zeros(2*N-1,1);
Ls(1:2:end) = li;
Ls(2:2:end) = lij;
xipm = [0,cumsum(Ls)']+spacing; % all boundary points
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
C = make_capacitance_finite(N,lij); % capacitance matrix
w_muller = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0,lij,xm,xp); % subwavelength resonant frequencies
w_res = w_muller(real(w_muller)>=0); % positive subwavelength resonant frequencies
w_op = w_res(1)+ 0.0002; % operating frequency
w0 = real(w_res(1))+0.02; % quasifrequency of incident wave
k_op = w_op/v0; % operating wave number outside of the resonator
k0 = w0/v0; % wave number of incident frequency

% Define relevant functions
uin_l = @(x,t,n) exp(sqrt(-1)*((k0).*x+w0.*t)).*(x<xm(1)).*(n==0); % left incident wave
vin_l = @(x,n) exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % n-th mode of left incident wave
dx_vin_l = @(x,n) sqrt(-1)*k0*exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % derivative of left incident wave
uin_r = @(x,t,n) exp(sqrt(-1)*(-(k0).*x+w0.*t)).*(x>xp(end)).*(n==0); % right incident wave
vin_r = @(x,n) exp(sqrt(-1)*(-(k0).*x)).*(x>xp(end)).*(n==0); % n-th mode of right incident wave
dx_vin_r = @(x,n) -sqrt(-1)*k0*exp(sqrt(-1)*(-(k0).*x)).*(x>xp(end)).*(n==0); % derivative of right incident wave
G = @(k,x) exp(sqrt(-1)*k*abs(x))./(2*sqrt(-1)*k); % Green's function

% Define evaluation points
len_xs = 800;
len_zs = 80;
xs = linspace(xm(1)-1,xp(end)+1,len_xs);
zs = zeros(N,len_zs);
for i = 1:N
    zs(i,:) = linspace(xm(i),xp(i),len_zs);
end

%% Compute the scattered wave field

% Compute solution coefficients
MatcalA = getMatcalA(N,lij,xm,xp,k_tr,w_op,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
MatcalF = getF_lr(k_tr, N, delta, xm, xp, dx_vin_l, dx_vin_r); % vector \mathcal{F}
sol = MatcalA\MatcalF; % solve for the interior coefficients, vector \mathbf{w}

us_eval_x = zeros(1,len_xs);
us_eval_z = zeros(1,len_zs);
usx = @(x) N*(uin_l(x,t,0)+uin_r(x,t,0)) + get_us_lr(x, t, N, xm, xp, lij, k_tr, v0, w_op, Omega, rs, ks, vr, sol, w_res, k0, vin_l, vin_r); % scattered wave field as a function of x for fixed time t, according to formula (31)

for i = 1:len_xs
    us_eval_x(i) = usx(xs(i));
end

% Create plot
fig = figure();
fig.Position = [263,725,982,352];
subplot(1,2,1)
set(gca,'FontSize',14)
hold on
plot(xs,real(us_eval_x),'.','DisplayName','Exact',markersize=8,linewidth=2)
xlim([xs(1) xs(end)])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel(strcat('Re$(u^{\mathrm{sc}}(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)

subplot(1,2,2)
set(gca,'FontSize',14)
hold on
plot(xs,imag(us_eval_x),'.','DisplayName','Exact',markersize=8,linewidth=2)
xlim([xs(1) xs(end)])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel(strcat('Im$(u^{\mathrm{sc}}(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)

for i = 1:N
    for j = 1:len_zs
        us_eval_z(j) = usx(zs(i,j));
    end
    subplot(1,2,1)
    plot(zs(i,:),real(us_eval_z),'r.','DisplayName','Exact',markersize=8,linewidth=2)
    subplot(1,2,2)
    plot(zs(i,:),imag(us_eval_z),'r.','DisplayName','Exact',markersize=8,linewidth=2)
end

