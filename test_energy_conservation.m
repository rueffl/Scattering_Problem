clear
format long

% Settings for the material's structure
k_tr = 4; % truncation parameters as in remark 3.3
N = 1; % number of the resonator
spacing = 10; lij = ones(1,N-1).*spacing; % spacing between the resonators
len = 1; li = ones(1,N).*len; % length of the resonator
L = sum(li)+sum(lij); % length of the unit cell
Ls = zeros(2*N-1,1);
Ls(1:2:end) = li;
Ls(2:2:end) = lij;
xipm = [0,cumsum(Ls)']-(len*(N/2)+spacing*(N-1)/2); % all boundary points, make sure the resonators are aligned symmetrically wrt 0
xm = xipm(1:2:end); % LHS boundary points
xp = xipm(2:2:end); % RHS boundary points
delta = 0.00001; % small contrast parameter
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
if N == 1
    w_res = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0); % non-zero subwavelength resonant frequency
else
    C = make_capacitance_finite(N,lij); % capacitance matrix
    w_muller = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0,lij,xm,xp); % subwavelength resonant frequencies
    w_res = w_muller(real(w_muller)>=0); % positive subwavelength resonant frequencies
end
w_op = w_res(1)+ 0.0002; % operating frequency
k_op = w_op/v0; % operating wave number outside of the resonator
w0 = w_op; k0 = k_op;
G = @(k,x) exp(sqrt(-1)*k*abs(x))./(2*sqrt(-1)*k); % Green's function
dx_G = @(k,x) x*exp(sqrt(-1)*k*abs(x))./(2*sqrt(-1)*abs(x)); % Derivative of Green's function

% fig = figure();
% hold on
% c_map = parula(2*k_tr+1+2); % colors for plot


w_op = w0;
k_op = w_op/v0;
k0 = w0/v0; % wave number of incident frequency

% Define incident wave field

vin_l = @(x,n) 1*exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % n-th mode of left incident wave
dx_vin_l = @(x,n) 1*sqrt(-1)*k0*exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % derivative of left incident wave
% uin_l = @(x,t,n) 0; vin_l = @(x,n) 0; dx_vin_l = @(x,n) 0;
% uin_r = @(x,t,n) exp(sqrt(-1)*(-(k0).*x+w0.*t)).*(x>xp(end)).*(n==0); % right incident wave
% vin_r = @(x,n) exp(sqrt(-1)*(-(k0).*x)).*(x>xp(end)).*(n==0); % n-th mode of right incident wave
% dx_vin_r = @(x,n) -sqrt(-1)*k0*exp(sqrt(-1)*(-(k0).*x)).*(x>xp(end)).*(n==0); % derivative of right incident wave
uin_r = @(x,t,n) 0; vin_r = @(x,n) 0; dx_vin_r = @(x,n) 0;

% Build the linear system and solve it for the interior coefficients
MatcalA = getMatcalA(N,lij,xm,xp,k_tr,w_op,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
MatcalF = getF_lr(k_tr, N, delta, xm, xp, dx_vin_l, dx_vin_r); % vector \mathcal{F}
N_vin = getN_vin(k_tr, N, delta, xm, xp, w_op, Omega, v0, lij, vin_l, vin_r); % matrix-vector product \mathcal{N}v^{in}
RHS = MatcalF + N_vin;
sol = MatcalA\RHS; % solve for the interior coefficients, vector \mathbf{w}

% Compute the transmission and reflection coefficients
[tns,rns] = get_tr(k_tr, w_op, Omega, rs, ks, vr, xm, xp, sol, lij, vin_l, vin_r);
O = abs(tns).^2+abs(rns).^2;

% Create plot of the transmission and reflection coefficients
%     for i = 1:N
%         figure(i)
%         hold on
%     end





% xlabel('$\omega_{\mathrm{in}}$',Interpreter='latex',FontSize=18)
% ylabel('$\Lambda_n$',Interpreter='latex',FontSize=18)
% legend('show',interpreter='latex',fontsize=18)


