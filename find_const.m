% main script
format long

%% Set parameters

% Settings for the structure
k_tr = 4; % truncation parameters as in remark 3.3
N = 6; % number of the resonator
li = ones(1,N).*0.025; % length of the resonators
lij = ones(1,N-1).*500; % distance between the resonators
L = sum(li)+sum(lij); % length of the unit cell
Ls = zeros(2*N-1,1);
Ls(1:2:end) = li;
Ls(2:2:end) = lij;
xipm = [0,cumsum(Ls)'];
xm = xipm(1:2:end);
xp = xipm(2:2:end);
z = (xm+xp)./2;
delta = 0.0001; % small contrast parameter
vr = 1;
v0 = 1;

% Settings for modulation
Omega = 0.03; % modulation frequency
T = 2*pi/Omega;
phase_kappa = zeros(1,N); % modulation phases of kappa
phase_rho = zeros(1,N); % modulation phases of rho
for i = 1:(N-1)
    phase_kappa(i+1) = pi/i;
    phase_rho(i+1) = pi/i;
end
epsilon_kappa = 0.2; % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho
rs = [];
ks = [];
for j = 1:N
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
    ks = [ks; ks_j];
    rs = [rs; rs_j];
end

%green's function
G = @(x,kn) exp(sqrt(-1)*kn*abs(x))./(2*sqrt(-1)*kn);

if N > 1
    C = make_capacitance_finite(N,lij);
    w_out = get_capacitance_approx(epsilon_kappa,epsilon_rho,li,Omega,phase_rho,phase_kappa,delta,C);
    w_muller = w_out(real(w_out)>=0);
else
    w_muller = [0];
end

w = w_muller(1) + 0.0001;
k = w/v0;
x = L+100;

gamma_s = zeros(2*k_tr+1,N);
const = zeros(2*k_tr+1,N);
for n = -k_tr:k_tr
    for j = 1:N

        w0 = mean(w_muller);
        k_0 = w0/v0;

        A = getMatcalA(N, lij, xm, xp, k_tr, w0, Omega, rs, ks, vr, delta, v0);
        F = getF(k_tr, N, delta, k, k_0, xm);
        sol = linsolve(A,F);

        gamma_s(n+k_tr+1,j) = operator_S(x, N, xm, xp, lij, k_tr, k, w, Omega, rs, ks, vr, sol, n)*(w-w0);
        const(n+k_tr+1,j) = gamma_s(n+k_tr+1,j)./(G(x-z(j),(w+n*Omega)./v0).*exp(sqrt(-1).*k.*z(j)));

    end
end










