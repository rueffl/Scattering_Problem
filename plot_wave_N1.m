% Plot the scattered wave in 1D for the single-resonator case

format long

% Settings for the structure
k_tr = 2; % truncation parameters as in remark 3.3
N = 1; % number of the resonator
spacing = 1; lij = ones(1,N-1).* spacing; % spacing between the resonators
len = 1; li = ones(1,N).*len; % length of the resonator
xm = 0; xp = xm+li(1);
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
epsilon_kappa = 0; % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho
rs = [];
ks = [];
for j = 1:N
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
    ks = [ks; ks_j];
    rs = [rs; rs_j];
end
% resonant frequencies
w_out = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0);
w_res = zeros(1,2); w_res(1,2) = w_out; 
w0 = w_out + 0.0002; % operating frequency
w = w0 + 0.00000001; % quasifrequency of incident wave

% Compute scattered wave
ls = 100; % number of evaluation points
xs = linspace(xm(1)-1, xp(N)+1,ls); % evaluation points
us = zeros(ls,1);

t = 0; % time
k = w/v0; % wave number of incident wave
k_0 = w0/v0; % wave number of operating wave

for i = 1:ls % iterate over evaluation points
    for j = 1:N % iterate over resonaros

        A = getMatcalA(N, lij, xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0); % matrix \mathcal{A}
        F = getF(k_tr, N, delta, k, k_0, xm); % vector \mathcal{F}
        sol = linsolve(A,F); % interior coefficients
        us(i) = us(i) + get_us(xs(i), t, N, xm, xp, lij, k_tr, v0, w, Omega, rs, ks, vr, sol, w0); % scattered wave at xs(i)

    end
end

% create plot
fig1 = figure()
hold on
plot(xs,real(us),'-',LineWidth=5)
title('Real Part of the Scattered Wave',fontsize=40)

fig2 = figure()
hold on
plot(xs,imag(us),'-',LineWidth=5)
title('Imaginary Part of the Scattered Wave',fontsize=40)

