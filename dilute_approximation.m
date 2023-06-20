% main script
format long

%% One Resonator
% Settings for the structure
k_tr = 3; % truncation parameters as in remark 3.3
N = 1; % number of the resonator
li = ones(1,N).*0.01; % length of the resonators
lij = ones(1,N-1).*100; %lij(2:2:end) = lij(2:2:end)+1; % distance between the resonators
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
epsilon_kappa = 0; % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho
rs = [];
ks = [];
for j = 1:N
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j)),1,epsilon_rho*exp(-1i*phase_rho(j))];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j)),1,epsilon_kappa*exp(-1i*phase_kappa(j))];
    ks = [ks; ks_j];
    rs = [rs; rs_j];
end

ls = 100;
us = zeros(ls,1);
xs = linspace(xm(1)-1, xp(N)+1,ls);

% Find quasifrequencies
w_muller = muller(0,N,lij,xm,xp,k_tr,Omega,rs,ks,vr,delta,v0);

t = 0;
w = w_muller(1) + 0.0001;
k = w/v0;
w0 = w_muller;
k_0 = w0/v0;
A = getMatcalA(N, lij, xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0);
F = getF(k_tr, N, delta, k, k_0, xm);
sol = linsolve(A,F);

all_cn = zeros(1,2*k_tr+1);
for n = -k_tr:k_tr
    kn = (w+n*Omega)/v0;
    lambdas = zeros(1,2*k_tr+1);
    fn = zeros(1,2*k_tr+1);
    for j = -k_tr:k_tr
        C = getC(k_tr, 1, w, Omega, rs, ks, vr);
        [fi,lls] = eig(C,'vector');
        lls = sqrt(lls);
        lambdas(j+k_tr+1) = lls(k_tr-j+1);
        if k_tr-n+1 < 1 || k_tr-n+1 > 2*k_tr+1
            fn(j+k_tr+1) = 0;
        else
            fn(j+k_tr+1) = fi(k_tr-n+1,k_tr-j+1);
        end
    end
    
    all_cn(n+k_tr+1) = get_cn(sol,kn,fn,lambdas,k_tr,z,1);
end

%green's function
G = @(x,kn) exp(sqrt(-1)*kn*abs(x))./(2*sqrt(-1)*kn);

%plot scattered wave
xs = linspace(z-2,z+2,400);
us_x = zeros(1,length(xs));
j = 1;
figure()
hold on
for x = xs
    us = 0;
    for n = -k_tr:k_tr
        us = us + all_cn(n+k_tr+1).*G(x,(w+n*Omega)./v0).*exp(sqrt(-1)*(w0+n*Omega).*t)./(w-w0);
    end
   us_x(j) = us;
   j = j + 1;
end
plot(xs,us_x,'b-')

%compare with exact formula for small resonator
[us, xs] = prepare_plot(k_tr,li,lij,L,N,Omega,delta,phase_kappa,phase_rho,epsilon_kappa,epsilon_rho,xm,xp,vr,v0);


%% N Resonators
% Settings for the structure
k_tr = 3; % truncation parameters as in remark 3.3
N = 6; % number of the resonator
li = ones(1,N).*0.001; % length of the resonators
lij = ones(1,N-1).*100; %lij(2:2:end) = lij(2:2:end)+1; % distance between the resonators
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
epsilon_kappa = 0; % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho
rs = [];
ks = [];
for j = 1:N
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j)),1,epsilon_rho*exp(-1i*phase_rho(j))];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j)),1,epsilon_kappa*exp(-1i*phase_kappa(j))];
    ks = [ks; ks_j];
    rs = [rs; rs_j];
end

% Find quasifrequencies
w_muller = zeros(N,1);

% Compute static case
C = make_capacitance_finite(N,lij);
w_static = get_capacitance_approx(0,0,li,Omega,phase_rho,phase_kappa,delta,C);

% Compute with muller's method 
for i = 1:2*N
    initial_guess = w_static(i);
    w_muller(i) = muller(initial_guess,N,lij,xm,xp,k_tr,Omega,rs,ks,vr,delta,v0);
end

t = 0;
w = w_muller(1) + 0.0001;
k = w/v0;
w0 = w_muller;
k_0 = w0/v0;
A = getMatcalA(N, lij, xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0);
F = getF(k_tr, N, delta, k, k_0, xm);
sol = linsolve(A,F);

% Define coefficients alpha_n
alphas = zeros(1,2*k_tr+1);
for n = -k_tr:k_tr
    kn = (w+n*Omega)/v0;
    lambdas = zeros(1,2*k_tr+1);
    fn = zeros(1,2*k_tr+1);
    for j = -k_tr:k_tr
        C = getC(k_tr, 1, w, Omega, rs, ks, vr);
        [fi,lls] = eig(C,'vector');
        lls = sqrt(lls);
        lambdas(j+k_tr+1) = lls(k_tr-j+1);
        if k_tr-n+1 < 1 || k_tr-n+1 > 2*k_tr+1
            fn(j+k_tr+1) = 0;
        else
            fn(j+k_tr+1) = fi(k_tr-n+1,k_tr-j+1);
        end
    end
    
    alphas(n+k_tr+1) = get_cn(sol,kn,fn,lambdas,k_tr,z,1)*exp(-sqrt(-1)*k_0);
end


