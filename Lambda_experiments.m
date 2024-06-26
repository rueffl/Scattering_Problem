%% Set structure setting
clear all
format long

% Settings for the material's structure
k_tr = 4; % truncation parameters as in remark 3.3
k_tr_n = k_tr;
k_tr_m = 0;
N = 1; % number of the resonator
spacing = 10; lij = ones(1,N-1).*spacing; % spacing between the resonators
len = 0.5; li = ones(1,N).*len; % length of the resonator
L = sum(li)+sum(lij); % length of the unit cell
Ls = zeros(2*N-1,1);
Ls(1:2:end) = li;
Ls(2:2:end) = lij;
xipm = [0,cumsum(Ls)']+spacing; % all boundary points
xm = xipm(1:2:end); % LHS boundary points
xp = xipm(2:2:end); % RHS boundary points
z = (xm+xp)./2; % centers of resonators
delta = 0.0001; % small contrast parameter
t = 0; % time

vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators


%% Compute the frequency-scattering coefficients for the static case

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
w_res = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0); % non-zero subwavelength resonant frequency
w_op = w_res(1)+ 0.0002; % operating frequency
w_in = w_res(1); % quasifrequency of incident wave
k_op = w_op/v0; % operating wave number outside of the resonator
k_in = w_in/v0; % wave number of incident frequency

% Solve the linear system for the solution coefficients
MatcalA = getMatcalA(N,lij,xm,xp,k_tr,w_op,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
MatcalF = getF(k_tr, N, delta, k_op, k_in, xm); % vector \mathcal{F}
sol = MatcalA\MatcalF; % solve for the interior coefficients, vector \mathbf{w}

x = 10;
Lambdas = get_Lambdas_N1(z,xm,xp,k_tr,w_op,Omega,rs,ks,vr,v0,delta,k_op,k_in,x,1,sol);

% Create plot
figure()
hold on
plot(-k_tr:k_tr,abs(Lambdas),'*')

%% Iterate over modulation amplitude and plot the frequency-scattering coefficient

% Setting for the material's time-modulation
Omega = 0.03; % modulation frequency
T = 2*pi/Omega;
phase_kappa = zeros(1,N); % modulation phases of kappa
phase_rho = zeros(1,N); % modulation phases of rho
for i = 1:(N-1)
    phase_kappa(i+1) = pi/i;
    phase_rho(i+1) = pi/i;
end
all_epsk = linspace(0,0.9,10); % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho

i = 1;
Lambdas1 = zeros(length(all_epsk),2*k_tr+1);
Lambdas2 = zeros(length(all_epsk),2*k_tr+1);
c_map = parula(length(all_epsk)+2);
figure();
hold on
for epsilon_kappa = all_epsk

    rs = []; % Fourier coefficients of 1/rho
    ks = []; % Fourier coefficients of 1/kappa
    for j = 1:N
        rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end

    % Calculate subwavelength resonant frequency
    w_res = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0); % non-zero subwavelength resonant frequency
    w_op = w_res(1)+ 0.0002; % operating frequency
    w_in = w_res(1); % quasifrequency of incident wave
    k_op = w_op/v0; % operating wave number outside of the resonator
    k_in = w_in/v0; % wave number of incident frequency

    MatcalA = getMatcalA(N,lij,xm,xp,k_tr,w_op,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
    MatcalF = getF(k_tr, N, delta, k_op, k_in, xm); % vector \mathcal{F}
    sol = MatcalA\MatcalF; % solve for the interior coefficients, vector \mathbf{w}

    x1 = 1;
    Lambdas1(i,:) = get_Lambdas_N1(z,xm,xp,k_tr,w_op,Omega,rs,ks,vr,v0,delta,k_op,k_in,x1,sol);
    x2 = 100;
    Lambdas2(i,:) = get_Lambdas_N1(z,xm,xp,k_tr,w_op,Omega,rs,ks,vr,v0,delta,k_op,k_in,x2,sol);

    subplot(1,2,1)
    hold on
    plot(-k_tr:k_tr,abs(Lambdas1(i,:)),'*','Color',c_map(i,:),'MarkerSize',10,'LineWidth',2,'DisplayName',strcat('$\varepsilon_{\kappa}=$',num2str(epsilon_kappa)))
    subplot(1,2,2)
    hold on
    plot(-k_tr:k_tr,abs(Lambdas2(i,:)),'*','Color',c_map(i,:),'MarkerSize',10,'LineWidth',2,'DisplayName',strcat('$\varepsilon_{\kappa}=$',num2str(epsilon_kappa)))
    i = i+1;
end

subplot(1,2,1)
ylabel(strcat('$|\Lambda_n($',num2str(x1),'$,$',num2str((xm+xp)/2),'$)|$'),'Interpreter','latex')
xlabel('$n$','Interpreter','latex')
subplot(1,2,2)
ylabel(strcat('$|\Lambda_n($',num2str(x2),'$,$',num2str((xm+xp)/2),'$)|$'),'Interpreter','latex')
xlabel('$n$','Interpreter','latex')
legend('show',interpreter='latex',fontsize=18,location='southoutside')


%% Plot \Lambda_n as a function of \ell

all_len = linspace(0.01,1,100);
k_tr = 2;

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
t = 0; % time

all_Ln1 = zeros(2*k_tr+1,length(all_len)); 
all_Ln2 = zeros(2*k_tr+1,length(all_len)); 
j = 1;
for len = all_len

    li = ones(1,N).*len; % length of the resonator
    xm = 10; xp = xm+li(1);

    % resonant frequencies
    w_out = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0);
    w_res = zeros(1,2); w_res(1,2) = w_out; 
    w_op = w_out + 0.0002; % operating frequency
    w_in = real(w_op) + 0.000001; % quasifrequency of incident wave
    k_in = w_in/v0; % wave number of incident wave
    k_op = w_op/v0; % wave number of operating wave

    MatcalA = getMatcalA(N,lij,xm,xp,k_tr,w_op,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
    MatcalF = getF(k_tr, N, delta, k_op, k_in, xm); % vector \mathcal{F}
    sol = MatcalA\MatcalF; % solve for the interior coefficients, vector \mathbf{w}

    x1 = 1;
    all_Ln1(:,j) = get_Lambdas_N1(xm,xp,k_tr,w_op,Omega,rs,ks,vr,v0,delta,k_op,k_in,x1,sol);
    x2 = 100;
    all_Ln2(:,j) = get_Lambdas_N1(xm,xp,k_tr,w_op,Omega,rs,ks,vr,v0,delta,k_op,k_in,x2,sol);
    j = j+1; 

end

% create plot of absolute value
c_map = parula(2*(2*k_tr+1)+1); iR = 1;
fig = figure();
fig.Position = [996,561,611,401];
for n = 1:(2*k_tr+1)
    loglog(all_len,abs(all_Ln1(n,:)),'-','Color',c_map(iR,:),'DisplayName',strcat('$n=$',num2str(n-(k_tr+1))),markersize=8,linewidth=2)
    hold on 
    loglog(all_len,abs(all_Ln2(n,:)),'--','Color',c_map(iR,:),'HandleVisibility','off',markersize=8,linewidth=2) 
    iR = iR+2;
end
% add legends etc.
legend('show',interpreter='latex',fontsize=18,location='southoutside')
xlabel('$\ell$',fontsize=18,Interpreter='latex')
ylabel(strcat('$|\Lambda_n(x,$',num2str((xm+xp)/2),'$)|$'),'Interpreter','latex','FontSize',18)

