% main script
format long

%% Set parameters

% Settings for the structure
N = 6; % number of the resonator
li = [1,1,1,1,1,1]; % length of the resonators
lij = [1,2,1,2,1]; % distance between the resonators
L = sum(li)+sum(lij); % length of the unit cell
Ls = zeros(2*N-1,1);
Ls(1:2:end) = li;
Ls(2:2:end) = lij;
xipm = [0,cumsum(Ls)'];
xm = xipm(1:2:end);
xp = xipm(2:2:end);
% xm = [0,li+lij(1:end-1)]; % define the boundary points x_minus and x_plus
% xp = xm + li; 
delta = 0.000001; % small contrast parameter

% Settings for modulation
Omega = 0.02; % modulation frequency
T = 2*pi/Omega;
epsilon_kappa = 0.1; % modulation amplitudes
epsilon_rho = 0;
phase_kappa = [0,pi/2,pi,pi/4,pi/6,pi/8]; % modulation phases
phase_rho = [0,pi/2,pi,pi/4,pi/6,pi/8];

k_tr = 4; % truncation parameters as in remark 3.3
vr = 1;
v0 = 1;

% Fourier coefficients of rhos and kappas
rs = [];
ks = [];
for j = 1:N
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j)),1,epsilon_rho*exp(-1i*phase_rho(j))];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j)),1,epsilon_kappa*exp(-1i*phase_kappa(j))];
    ks = [ks; ks_j];
    rs = [rs; rs_j];
end

%% Find quasifrequencies
w_muller = zeros(N,1);

% compute static case
C = make_capacitance_finite(N,lij);
w_static= get_capacitance_approx(0,0,li,Omega,phase_rho,phase_kappa,delta,C);

% % plot eigenvalues
% ev = [];
% wspan = linspace(-0.01,0.01,1000);
% for w = wspan
%    MatcalA = getMatcalA(N,lij,xm,xp,k_tr,w,Omega,rs,ks,vr,delta,v0);
%    ev = [ev svds(MatcalA,1,"smallest")]; %minev(MatcalA)];
% end
% figure
% plot(wspan,ev);

% compute with muller's method 
for i = 1:N
    initial_guess = w_static(i);
%     w_muller(i) = muller(initial_guess,N,lij,L,xm,xp,k_tr,Omega,rs,ks,vr,delta,v0);
    [w_muller(i),p0] = search_vector_finite(initial_guess,N,lij,xm,xp,k_tr,Omega,rs,ks,vr,delta,v0);
end


%% Compute scattered wave

ls = 100;
us = zeros(ls,1);
xs = linspace(xm(1)-1, xp(N)+1,ls);
% xs = linspace(xm(1)-1, xp(1),ls);

t = 0;

for i = 1:ls
    for j = 1:N

        w0 = w_muller(j);
        k_0 = w0/v0;
        w = w0 + 0.0001;
        k = w/v0;

        A = getMatcalA(N, lij, xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0);
        F = getF(k_tr, N, delta, k, k_0, xm);
        sol = linsolve(A,F);
        us(i) = us(i) + get_us(xs(i), t, N, xm, xp, lij, k_tr, k, w, Omega, rs, ks, vr, sol, w0);

    end
end
fig1 = figure()
plot(xs,real(us),'-')
title('Real Part of the Scattered Wave')
xlabel('$x$',Interpreter='latex')
ylabel('Re$\left(u^s\right)$',Interpreter='latex')

fig2 = figure()
plot(xs,imag(us),'-')
title('Imaginary Part of the Scattered Wave')
xlabel('$x$',Interpreter='latex')
ylabel('Im$\left(u^s\right)$',Interpreter='latex')




