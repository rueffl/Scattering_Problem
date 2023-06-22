format long

% Settings for the structure
k_tr = 4; % truncation parameters as in remark 3.3
N = 6; % number of the resonator
li = ones(1,N).*0.025; % length of the resonators
lij_s = linspace(0.5,500000,100); % distance between the resonators
delta = 0.0001; % small contrast parameter


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


all_w = zeros(length(lij_s),2*N);
j = 1;
figure()
hold on
for len = lij_s
    lij = ones(1,N-1).*len;
    C = make_capacitance_finite(N,lij);
    w_out = get_capacitance_approx(epsilon_kappa,epsilon_rho,li,Omega,phase_rho,phase_kappa,delta,C);
    plot(ones(1,2*N).*len,w_out,'k.',markersize=8)
    all_w(j,:) = w_out;
    j = j+1;
end
xlabel('$\ell_{i(i+1)}$',interpreter='latex',fontsize=20)
ylabel('$\omega_j$',interpreter='latex',fontsize=20)

