% main script
format long

%% Set parameters

% Settings for the structure
k_tr = 4; % truncation parameters as in remark 3.3
N = 6; % number of the resonator
li = ones(1,N); % length of the resonators
lij = ones(1,N-1); lij(2:2:end) = lij(2:2:end)+1; % distance between the resonators
L = sum(li)+sum(lij); % length of the unit cell
Ls = zeros(2*N-1,1);
Ls(1:2:end) = li;
Ls(2:2:end) = lij;
xipm = [0,cumsum(Ls)'];
xm = xipm(1:2:end);
xp = xipm(2:2:end);
delta = 0.000001; % small contrast parameter
vr = 1;
v0 = 1;

% Settings for modulation
Omega = 0.02; % modulation frequency
T = 2*pi/Omega;
phase_kappa = zeros(1,N); % modulation phases of kappa
phase_rho = zeros(1,N); % modulation phases of rho
for i = 1:(N-1)
    phase_kappa(i+1) = pi/i;
    phase_rho(i+1) = pi/i;
end


epsilon_kappa = 0; % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho
[us00, xs] = prepare_plot(k_tr,li,lij,L,N,Omega,delta,phase_kappa,phase_rho,epsilon_kappa,epsilon_rho,xm,xp,vr,v0);

epsilon_kappa = 0.2; % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho
[us20, xs] = prepare_plot(k_tr,li,lij,L,N,Omega,delta,phase_kappa,phase_rho,epsilon_kappa,epsilon_rho,xm,xp,vr,v0);

epsilon_kappa = 0; % modulation amplitude of kappa
epsilon_rho = 0.2; % modulation amplitude of rho
[us02, xs] = prepare_plot(k_tr,li,lij,L,N,Omega,delta,phase_kappa,phase_rho,epsilon_kappa,epsilon_rho,xm,xp,vr,v0);

epsilon_kappa = 0.2; % modulation amplitude of kappa
epsilon_rho = 0.2; % modulation amplitude of rho
[us22, xs] = prepare_plot(k_tr,li,lij,L,N,Omega,delta,phase_kappa,phase_rho,epsilon_kappa,epsilon_rho,xm,xp,vr,v0);

fig1 = figure(1)
hold on
plot(xs,real(us00),'-')
plot(xs,real(us20),'-')
plot(xs,real(us02),'-')
plot(xs,real(us22),'-')
legend('$\varepsilon_{\kappa}=0,\,\varepsilon_{\rho}=0$','$\varepsilon_{\kappa}=0.2,\,\varepsilon_{\rho}=0$', ...
    '$\varepsilon_{\kappa}=0,\,\varepsilon_{\rho}=0.2$','$\varepsilon_{\kappa}=0.2,\,\varepsilon_{\rho}=0.2$', ...
    interpreter='latex')
title('Real Part of the Scattered Wave')

fig2 = figure(2)
hold on
plot(xs,imag(us00),'-')
plot(xs,imag(us20),'-')
plot(xs,imag(us02),'-')
plot(xs,imag(us22),'-')
legend('$\varepsilon_{\kappa}=0,\,\varepsilon_{\rho}=0$','$\varepsilon_{\kappa}=0.2,\,\varepsilon_{\rho}=0$', ...
    '$\varepsilon_{\kappa}=0,\,\varepsilon_{\rho}=0.2$','$\varepsilon_{\kappa}=0.2,\,\varepsilon_{\rho}=0.2$', ...
    interpreter='latex')
title('Imaginary Part of the Scattered Wave')

