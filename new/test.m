format long

%% Compute \Lambda_n as a function of \ell

% Settings for the structure
k_tr = 2; % truncation parameter
N = 1; % number of the resonator
lij = ones(1,N-1); % spacing between the resonators
all_len = linspace(0.01,1,100);
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
t = 0; % time

all_Ln = zeros(2*k_tr+1,length(all_len)); 
j = 1;
for len = all_len
    li = ones(1,N).*len; % length of the resonator
    xm = 0; xp = xm+li(1);
    % resonant frequencies
    w_out = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0);
    w_res = zeros(1,2); w_res(1,2) = w_out; 
    w0 = w_out + 0.0002; % operating frequency
    w = real(w0) + 0.000001; % quasifrequency of incident wave
    k = w/v0; % wave number of incident wave
    k_0 = w0/v0; % wave number of operating wave
    all_Ln(:,j) = get_Lambdas(xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0, k, k_0, (xm+xp)/2, [])';
    j = j+1; 
end

% create plot of absolute value
c_map = parula(2*(2*k_tr+1)+1); iR = 1;
fig = figure();
fig.Position = [996,561,611,401];
for n = 1:(2*k_tr+1)
    loglog(all_len,abs(all_Ln(n,:)),'-','Color',c_map(iR,:),'DisplayName',strcat('$n=$',num2str(n-(k_tr+1))),markersize=8,linewidth=2)
    hold on  
    iR = iR+2;
end
% add legends etc.
legend('show',interpreter='latex',fontsize=18,location='southoutside')
xlabel('$\ell$',fontsize=18,Interpreter='latex')
ylabel('$|\Lambda_n(\ell)|$',fontsize=18,Interpreter='latex')

%% Compute \Lambda_n as a function of n

% Settings for the structure
k_tr = 10; % truncation parameters as in remark 3.3
N = 1; % number of the resonator
lij = ones(1,N-1); % spacing between the resonators
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
all_epsk = linspace(0,0.9,10);
% all_epsk = [0];

all_Ln = zeros(2*k_tr+1,length(all_epsk));
c = 1;

for epsilon_kappa = all_epsk % iterate over modulation amplitudes of kappa

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
    w = real(w0) + 0.00000001; % quasifrequency of incident wave
    k = w/v0; % wave number of incident wave
    k_0 = w0/v0; % wave number of operating wave
    t = 0; % time
    
    all_Ln(:,c) = get_Lambdas(xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0, k, k_0, (xm+xp)/2, [])'; % \Lambda_n for all n=-K,...,K
    c = c+1;

end

% create plot of imaginary and real parts
c_map = parula(length(all_epsk)+1);
figure();
for i = 1:length(all_epsk)
    subplot(1,2,1)
    plot(-k_tr:k_tr,real(all_Ln(:,i)),'*','Color',c_map(i,:),'DisplayName',strcat('$\varepsilon_{\kappa}=$',num2str(all_epsk(i))),markersize=8,linewidth=2)
    hold on
    subplot(1,2,2)
    plot(-k_tr:k_tr,imag(all_Ln(:,i)),'*','Color',c_map(i,:),'DisplayName',strcat('$\varepsilon_{\kappa}=$',num2str(all_epsk(i))),markersize=8,linewidth=2)
    hold on    
end
% add legends etc.
subplot(1,2,1)
legend('show',interpreter='latex',fontsize=18,location='southoutside')
xticks([-k_tr:k_tr])
xlabel('$n$',fontsize=18,Interpreter='latex')
ylabel(strcat('Re$(\Lambda_n($',num2str(len),'$))$'),fontsize=18,Interpreter='latex')
subplot(1,2,2)
legend('show',interpreter='latex',fontsize=18,location='southoutside')
xticks([-k_tr:k_tr])
xlabel('$n$',fontsize=18,Interpreter='latex')
ylabel(strcat('Im$(\Lambda_n($',num2str(len),'$))$'),fontsize=18,Interpreter='latex')

% create plot of absolute value
figure();
for i = 1:length(all_epsk)
    plot(-k_tr:k_tr,abs(all_Ln(:,i)),'*','Color',c_map(i,:),'DisplayName',strcat('$\varepsilon_{\kappa}=$',num2str(all_epsk(i))),markersize=8,linewidth=2)
    hold on
end
% add legends etc.
legend('show',interpreter='latex',fontsize=18,location='southoutside')
xticks([-k_tr:k_tr])
xlabel('$n$',fontsize=18,Interpreter='latex')
ylabel(strcat('$|\Lambda_n($',num2str(len),'$)|$'),fontsize=18,Interpreter='latex')


%% Compute \Lambda_n as a function of \delta
% Settings for the structure
k_tr = 10; % truncation parameters as in remark 3.3
N = 1; % number of the resonator
lij = ones(1,N-1); % spacing between the resonators
len = 1; li = ones(1,N).*len; % length of the resonator
xm = 0; xp = xm+li(1);
all_deltas = [10^(-7),10^(-6),10^(-5),10^(-4),10^(-3),10^(-2),10^(-1)]; % small contrast parameter

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
t = 0; % time

all_Ln = zeros(2*k_tr+1,length(all_deltas)); 
j = 1;
for delta = all_deltas
    Omega = 4*sqrt(delta); % modulation frequency
    T = 2*pi/Omega;
    % resonant frequencies
    w_out = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0);
    w_res = zeros(1,2); w_res(1,2) = w_out; 
    w0 = w_out + 0.0002; % operating frequency
    w = real(w0) + 0.000001; % quasifrequency of incident wave
    k = w/v0; % wave number of incident wave
    k_0 = w0/v0; % wave number of operating wave
    all_Ln(:,j) = get_Lambdas(xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0, k, k_0)';
    j = j+1; 
end

% create plot of absolute value
figure();
for i = 1:length(all_deltas)
    semilogy(-k_tr:k_tr,abs(all_Ln(:,i)),'*','Color',c_map(i,:),'DisplayName',strcat('$\delta=$',num2str(all_deltas(i))),markersize=8,linewidth=2)
    hold on
end
% add legends etc.
legend('show',interpreter='latex',fontsize=18,location='southoutside')
xticks([-k_tr:k_tr])
xlabel('$n$',fontsize=18,Interpreter='latex')
ylabel(strcat('$|\Lambda_n($',num2str(len),'$)|$'),fontsize=18,Interpreter='latex')


%% Experiments for \Lambda_{mn} as a matrix

% Settings for the structure
k_tr = 2; % truncation parameter
N = 1; % number of the resonator
lij = ones(1,N-1); % spacing between the resonators
len = 0.5; li = ones(1,N).*len; % length of the resonators
xm = 0; xp = xm+li(1); % boundary points of the resonator
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
t = 0; % time

% resonant frequencies
w_out = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0);
w_res = zeros(1,2); w_res(1,2) = w_out; 
w0 = w_out + 0.0002; % operating frequency
w = real(w0) + 0.000001; % quasifrequency of incident wave
k = w/v0; % wave number of incident wave
k_0 = w0/v0; % wave number of operating wave
all_Lnm = get_Lambdas_mn(xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0, k_0);

% create plot of absolute value
pcolor(-k_tr:k_tr,-k_tr:k_tr,abs(all_Lnm))

















