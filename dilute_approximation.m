format long
%% 1 Resonator
% Settings for the structure
k_tr = 1; % truncation parameters as in remark 3.3
N = 1; % number of the resonator
li = ones(1,N).*0.0001; % length of the resonators
lij = ones(1,N-1).*5000; %lij(2:2:end) = lij(2:2:end)+1; % distance between the resonators
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
x = L+10;

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
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j))/2,1,epsilon_rho*exp(1i*phase_rho(j))/2];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))/2,1,epsilon_kappa*exp(1i*phase_kappa(j))/2];
    ks = [ks; ks_j];
    rs = [rs; rs_j];
end

% Find quasifrequencies
w_muller = zeros(N,1);

% Compute static case
guess = zeros(1,2);
guess(2) = -(v0*vr*log((delta^2*vr^2 + 2*delta*v0*vr + v0^2)/(v0 - delta*vr)^2)*sqrt(-1))/(li(1)*(v0 + vr));%-sqrt(-1)*vr*log(1+2*vr*delta/(v0-vr*delta));
w_res = zeros(1,2);

for j = 1:2*N
    w_res(j) = muller(guess(j),N,lij,xm,xp,k_tr,Omega,rs,ks,vr,delta,v0);
end
% w_res = w_res(2);

%green's function
G = @(x,kn) exp(sqrt(-1)*kn*abs(x))./(2*sqrt(-1)*kn);

t = 0;
w_0 = mean(w_res);
w = w_0 + 0.0001; %quasifrequency of incident wave
k = w/v0;
k_0 = w_0/v0; %wave number of incident wave
A = getMatcalA(N, lij, xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0);
F = getF(k_tr, N, delta, k, k_0, z);
sol = linsolve(A,F);

% Define coefficients alpha_n
alphas = zeros(N,2*k_tr+1);
all_cn = zeros(N,2*k_tr+1);
gamma_s = zeros(N,2*k_tr+1);
const = zeros(N,2*k_tr+1);
for res = 1:N
    C = getC(k_tr, res, w, Omega, rs, ks, vr);
    [fi,lambdas] = eig(C,'vector');
    lls = sqrt(lambdas);
    lambdas = lls;
%     lambdas = zeros(1,2*k_tr+1);
%     for j = -k_tr:k_tr
%         lambdas(j+k_tr+1) = lls(k_tr-j+1);
%     end
    for n = -k_tr:k_tr
        kn = (w+n*Omega)/v0;
        fn = zeros(1,2*k_tr+1);
        for j = -k_tr:k_tr
            if k_tr-n+1 < 1 || k_tr-n+1 > 2*k_tr+1
                fn(j+k_tr+1) = 0;
            else
                fn(j+k_tr+1) = fi(k_tr-n+1,k_tr-j+1);
            end
        end
        
        alphas(res,n+k_tr+1) = get_cn(sol,kn,fn,lambdas,k_tr,z(res),res)*exp(-sqrt(-1)*k_0*z(res));
        all_cn(res,n+k_tr+1) = get_cn(sol,kn,fn,lambdas,k_tr,z(res),res);

        gamma_s(res,n+k_tr+1) = operator_S(x, N, z, z, lij, k_tr, kn, w, Omega, rs, ks, vr, sol, n)*(w-w_0);
        const(res,n+k_tr+1) = gamma_s(res,n+k_tr+1)./(G(x-z(res),kn));%.*exp(sqrt(-1).*k_0.*z(res)));

    end
end

%incident wave
% u_in = @(x,t) exp(sqrt(-1).*(k.*x+))

ls = 200;
all_x = linspace(L,L+40,ls);
ux = zeros(ls,1);
for l = 1:ls
    x = all_x(l);
    for n = -k_tr:k_tr
        kn = (w+n*Omega)/v0;
        for j = 1:N
            ux(l) = ux(l) + all_cn(j,n+k_tr+1)*G(x-z(j),kn).*exp(sqrt(-1)*(w_0+n*Omega)*t)./(w-w_0);
        end
    end
end

%compare with exact formula
[us, xs] = prepare_plot(k_tr,li,lij,L,N,Omega,delta,phase_kappa,phase_rho,epsilon_kappa,epsilon_rho,xm,xp,vr,v0,all_x);
us = us./(w-w_0);

figure()
hold on 
plot(all_x,real(us),'r-')
plot(all_x,real(ux),'c--')
legend('Exact Formula','Dilute Formula')
figure()
hold on 
plot(all_x,imag(us),'r-')
plot(all_x,imag(ux),'c--')
legend('Exact Formula','Dilute Formula')

rel_err = abs(us-ux)./abs(us);
figure()
plot(all_x,rel_err,'.')
legend('Relative Error')
figure()
plot(all_x,abs(us-ux),'.')
legend('Absolute Error')

%% N Resonators
% Settings for the structure
k_tr = 1; % truncation parameters as in remark 3.3
N = 2; % number of the resonator
li = ones(1,N).*0.025; % length of the resonators
lij = ones(1,N-1).*5000; %lij(2:2:end) = lij(2:2:end)+1; % distance between the resonators
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
x = L+1;

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
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j))/2,1,epsilon_rho*exp(1i*phase_rho(j))/2];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))/2,1,epsilon_kappa*exp(1i*phase_kappa(j))/2];
    ks = [ks; ks_j];
    rs = [rs; rs_j];
end

% Compute quasifrequencies
C = make_capacitance_finite(N,lij);
w_out = get_capacitance_approx(epsilon_kappa,epsilon_rho,li,Omega,phase_rho,phase_kappa,delta,C);
w_res = w_out(real(w_out)>=0);

%green's function
G = @(x,kn) exp(sqrt(-1)*kn*abs(x))./(2*sqrt(-1)*kn);

t = 0;
w_0 = mean(w_res);
w = w_0 + 0.0001; %quasifrequency of incident wave
k = w/v0;
k_0 = w_0/v0; %wave number of incident wave
A = getMatcalA(N, lij, xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0);
F = getF(k_tr, N, delta, k, k_0, z);
sol = linsolve(A,F);

% Define coefficients alpha_n
alphas = zeros(N,2*k_tr+1);
all_cn = zeros(N,2*k_tr+1);
gamma_s = zeros(N,2*k_tr+1);
const = zeros(N,2*k_tr+1);
for res = 1:N
    C = getC(k_tr, res, w, Omega, rs, ks, vr);
    [fi,lambdas] = eig(C,'vector');
    lls = sqrt(lambdas);
    lambdas = lls;
%     lambdas = zeros(1,2*k_tr+1);
%     for j = -k_tr:k_tr
%         lambdas(j+k_tr+1) = lls(k_tr-j+1);
%     end
    for n = -k_tr:k_tr
        kn = (w+n*Omega)/v0;
        fn = zeros(1,2*k_tr+1);
        for j = -k_tr:k_tr
            if k_tr-n+1 < 1 || k_tr-n+1 > 2*k_tr+1
                fn(j+k_tr+1) = 0;
            else
                fn(j+k_tr+1) = fi(k_tr-n+1,k_tr-j+1);
            end
        end
        
        alphas(res,n+k_tr+1) = get_cn(sol,kn,fn,lambdas,k_tr,z(res),res)*exp(-sqrt(-1)*k_0*z(res));
        all_cn(res,n+k_tr+1) = get_cn(sol,kn,fn,lambdas,k_tr,z(res),res);

        gamma_s(res,n+k_tr+1) = operator_S(x, N, z, z, lij, k_tr, k, w, Omega, rs, ks, vr, sol, n)*(w-w_0);
        const(res,n+k_tr+1) = gamma_s(res,n+k_tr+1)./(G(x-z(res),(w+n*Omega)./v0));%.*exp(sqrt(-1).*k_0.*z(res)));

    end
end

%incident wave
% u_in = @(x,t) exp(sqrt(-1).*(k.*x+))

ux = 0;
for n = -k_tr:k_tr
    kn = (w+n*Omega)/v0;
    for j = 1:N
        ux = ux + all_cn(j,n+k_tr+1)*G(x-z(j),kn).*exp(sqrt(-1)*(w_0+n*Omega)*t)./(w-w_0);
    end
end

%compare with exact formula
[us, xs] = prepare_plot(k_tr,li,lij,L,N,Omega,delta,phase_kappa,phase_rho,epsilon_kappa,epsilon_rho,xm,xp,vr,v0,[x]);
us = us./(w-w_0);
