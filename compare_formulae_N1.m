%% Set structure setting
clear all
format long

% Settings for the material's structure
k_tr = 4; % truncation parameters as in remark 3.3
k_tr_n = k_tr;
k_tr_m = 0;
N = 1; % number of the resonator
spacing = 10; lij = ones(1,N-1).*spacing; % spacing between the resonators
len = 1; li = ones(1,N).*len; % length of the resonator
L = sum(li)+sum(lij); % length of the unit cell
Ls = zeros(2*N-1,1);
Ls(1:2:end) = li;
Ls(2:2:end) = lij;
xipm = [0,cumsum(Ls)']; % all boundary points
xm = xipm(1:2:end)+spacing; % LHS boundary points
xp = xipm(2:2:end)+spacing; % RHS boundary points
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
w_op = w_res+0.0002; % operating frequency
w0 = 0.00088; % quasifrequency of incident wave
k_op = w_op/v0; % operating wave number outside of the resonator
k0 = w0/v0; % wave number of incident frequency

% Define relevant functions
uin = @(x,t) exp(sqrt(-1)*(k0.*x+w0.*t)).*(x<xm(1)); % incident wave
G = @(k,x) exp(sqrt(-1)*k*abs(x))./(2*sqrt(-1)*k); % Green's function

% Define evaluation points
len_xs = 800;
len_zs = 80;
xs = linspace(xm(1)-1,xp(end)+1,len_xs);
zs = linspace(xm(1),xp(end),len_zs);

%% Compute the scattered wave using formula (31)

% Compute solution coefficients
MatcalA = getMatcalA(N,lij,xm,xp,k_tr,w_op,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
MatcalF = getF(k_tr, N, delta, k_op, k0, xm); % vector \mathcal{F}
sol = MatcalA\MatcalF; % solve for the interior coefficients, vector \mathbf{w}

us_eval_x = zeros(1,len_xs);
us_eval_z = zeros(1,len_zs);
usx = @(x) uin(x,t) + get_us(x, t, N, xm, xp, lij, k_tr, v0, w_op, Omega, rs, ks, vr, sol, w_res, k0); % scattered wave field as a function of x for fixed time t, according to formula (31)

for i = 1:len_xs
    us_eval_x(i) = usx(xs(i));
end
for i = 1:len_zs
    us_eval_z(i) = usx(zs(i));
end

% Create plot
fig = figure();
fig.Position = [263,725,982,352];
subplot(1,2,1)
set(gca,'FontSize',14)
hold on
plot(xs,real(us_eval_x),'-','DisplayName','Exact',markersize=8,linewidth=2)
plot(zs,real(us_eval_z),'r-','DisplayName','Exact',markersize=8,linewidth=2)
xlim([xs(1) xs(end)])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel(strcat('Re$(u^{\mathrm{sc}}(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)

subplot(1,2,2)
set(gca,'FontSize',14)
hold on
plot(xs,imag(us_eval_x),'-','DisplayName','Exact',markersize=8,linewidth=2)
plot(zs,imag(us_eval_z),'r-','DisplayName','Exact',markersize=8,linewidth=2)
xlim([xs(1) xs(end)])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel(strcat('Im$(u^{\mathrm{sc}}(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)


%% Compute the scattered wave using formula (68)

all_Ln1 = get_Lambdas(xm, xp, k_tr, w_op, Omega, rs, ks, vr, delta, v0, k_op, k0, z, []); % frequency scattering coefficients

us_eval_x = zeros(1,len_xs);
% usx = @(x) u_sc(x,t,all_Ln,w_op,Omega,v0,w0,k0,k_tr_n,k_tr_m,z,w_res); % scattered wave field as a function of x for fixed time t, according to formula (68)
usx = @(x) u_sc(x,t,rs,ks,vr,delta,k_op,k0,w_op,Omega,v0,w0,k_tr_n,k_tr_m,z,w_res);

for i = 1:len_xs
    us_eval_x(i) = usx(xs(i));
end
us_eval_z = usx(z);

% Create plot
fig2 = figure();
fig2.Position = [263,725,982,352];
subplot(1,2,1)
set(gca,'FontSize',14)
hold on
plot(xs,real(us_eval_x),'-','DisplayName','Exact',markersize=8,linewidth=2)
plot(z,real(us_eval_z),'r-','DisplayName','Exact',markersize=8,linewidth=2)
xlim([xs(1) xs(end)])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel(strcat('Re$(u^{\mathrm{sc}}(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)

subplot(1,2,2)
set(gca,'FontSize',14)
hold on
plot(xs,imag(us_eval_x),'-','DisplayName','Exact',markersize=8,linewidth=2)
plot(z,imag(us_eval_z),'r-','DisplayName','Exact',markersize=8,linewidth=2)
xlim([xs(1) xs(end)])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel(strcat('Im$(u^{\mathrm{sc}}(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)































