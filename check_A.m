%% N Resonators
% Settings for the structure
k_tr = 3; % truncation parameters as in remark 3.3
N = 4; % number of the resonator
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

p = @(w) svds(getMatcalA(N, lij, xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0),1,"smallest");
% ws = linspace(-Omega+10^(-4),Omega-10^(-4),200);
ws = linspace(-Omega/2,Omega/2,200);
ps = zeros(1,200);
for i = 1:200
    ps(i) = p(ws(i));
end
figure()
plot(ws,ps)
xlabel('$\omega$',Interpreter='latex',fontsize=14)
title(strcat('Smallest eigenvalue of $\mathcal{A}$ for $\varepsilon_{\kappa}= $ ',num2str(epsilon_kappa),', $\varepsilon_{\rho}= $ ',num2str(epsilon_rho),', $\Omega= $ ',num2str(Omega),', $K =$ ',num2str(k_tr)),Interpreter="latex")    

    
% compute static case
C = make_capacitance_finite(N,lij);
w_static(j) = get_capacitance_approx_spec(0,zeros(N,1),Omega,delta,li,C);
w_cap(j) = get_capacitance_approx_spec(epsilon_kappa,phase_kappa,Omega,delta,li,C);


hold on
for i = 1:sample_points
    plot(w_cap(i),4*10^(-7),'r*')
end
for i = 1:sample_points
    plot3(w_static(i),4*10^(-7),'go')
end