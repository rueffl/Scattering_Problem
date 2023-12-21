% clear
% format long

% Settings for the material's structure
k_tr = 4; % truncation parameters as in remark 3.3
N = 6; % number of the resonator
spacing = 10; lij = ones(1,N-1).*spacing; % spacing between the resonators
len = 2; li = ones(1,N).*len; % length of the resonator
L = sum(li)+sum(lij); % length of the unit cell
Ls = zeros(2*N-1,1);
Ls(1:2:end) = li;
Ls(2:2:end) = lij;
xipm = [0,cumsum(Ls)']-(len*(N/2)+spacing*(N-1)/2); % all boundary points, make sure the resonators are aligned symmetrically wrt 0
xm = xipm(1:2:end); % LHS boundary points
xp = xipm(2:2:end); % RHS boundary points
delta = 0.0001; % small contrast parameter
t = 0; % time

vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators

% Setting for the material's time-modulation
Omega = 0.03; % modulation frequency
T = 2*pi/Omega;
phase_kappa = zeros(1,N); % modulation phases of kappa
phase_rho = zeros(1,N); % modulation phases of rho
for i = 1:(N-1)
    phase_kappa(i+1) = pi/i;
    phase_rho(i+1) = pi/i;
end
epsilon_kappa = 0.3; % modulation amplitude of kappa
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
if N == 1
    w_res = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0); % non-zero subwavelength resonant frequency
else
    C = make_capacitance_finite(N,lij); % capacitance matrix
    w_muller = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0); % subwavelength resonant frequencies
    w_res = w_muller(real(w_muller)>=0); % positive subwavelength resonant frequencies
end

% %% Iterate over omega
% 
% all_w = linspace(delta,sqrt(delta),10);
% O = zeros(2*k_tr+1,length(all_w)); E = zeros(1,length(all_w)); i = 1;
% % fig = figure();
% % hold on
% c_map = parula(length(all_w)+2);
% 
% for w_op = all_w(1:end)
%     
%     k_op = w_op/v0; % operating wave number outside of the resonator
%     w0 = w_op; k0 = k_op;
%     G = @(k,x) exp(sqrt(-1)*k*abs(x))./(2*sqrt(-1)*k); % Green's function
%     dx_G = @(k,x) x*exp(sqrt(-1)*k*abs(x))./(2*sqrt(-1)*abs(x)); % Derivative of Green's function
%     
%     % Define incident wave field  
%     vin_l = @(x,n) 1*exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % n-th mode of left incident wave
%     dx_vin_l = @(x,n) 1*sqrt(-1)*k0*exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % derivative of left incident wave
%     % uin_l = @(x,t,n) 0; vin_l = @(x,n) 0; dx_vin_l = @(x,n) 0;
%     % uin_r = @(x,t,n) exp(sqrt(-1)*(-(k0).*x+w0.*t)).*(x>xp(end)).*(n==0); % right incident wave
%     % vin_r = @(x,n) exp(sqrt(-1)*(-(k0).*x)).*(x>xp(end)).*(n==0); % n-th mode of right incident wave
%     % dx_vin_r = @(x,n) -sqrt(-1)*k0*exp(sqrt(-1)*(-(k0).*x)).*(x>xp(end)).*(n==0); % derivative of right incident wave
%     uin_r = @(x,t,n) 0; vin_r = @(x,n) 0; dx_vin_r = @(x,n) 0;
%     
%     % Build the linear system and solve it for the interior coefficients
%     MatcalA = getMatcalA(N,lij,xm,xp,k_tr,w_op,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
%     MatcalF = getF_lr(k_tr, N, delta, xm, xp, dx_vin_l, dx_vin_r); % vector \mathcal{F}
%     N_vin = getN_vin(k_tr, N, delta, xm, xp, w_op, Omega, v0, lij, vin_l, vin_r); % matrix-vector product \mathcal{N}v^{in}
%     RHS = MatcalF + N_vin;
%     sol = MatcalA\RHS; % solve for the interior coefficients, vector \mathbf{w}
%     
%     % Compute the transmission and reflection coefficients
%     [tns,rns] = get_tr(k_tr, w_op, Omega, rs, ks, vr, xm, xp, sol, lij, vin_l, vin_r);
%     O(:,i) = abs(tns(:,1)).^2+abs(rns(:,end)).^2;
%     E(i) = get_energy(tns(:,1),rns(:,end),w_op,Omega,v0,xm(1),xp(end),k_tr,t);
%     i = i+1;
% 
%     % Create plot
%     semilogy(-k_tr:k_tr,O(:,i-1),'o','Color',c_map(i-1,:),'DisplayName',strcat('$\omega=$',num2str(w_op)),'MarkerSize',8,'LineWidth',2)
%     hold on
% 
% end
% 
% xlabel('$n$',Interpreter='latex',FontSize=18)
% ylabel('$|R_n|^2+|T_n|^2$',Interpreter='latex',FontSize=18)
% legend('show',interpreter='latex',fontsize=18)
% 
% figure()
% hold on
% plot(all_w,E(1:end),'o','MarkerSize',8,'LineWidth',2)
% for w = w_res
%     plot(w.*ones(1,100),linspace(min(E(1:end)),max(E(1:end)),100),'--r','MarkerSize',8,'LineWidth',2)
% end
% xlabel('$\omega$',Interpreter='latex',FontSize=18)
% ylabel('Total Energy',Interpreter='latex',FontSize=18)
% 
% %% Iterate over epsilon_kappa
% 
% all_epsk = linspace(0,0.99,11);
% O = zeros(2*k_tr+1,length(all_epsk)); E = zeros(1,length(all_epsk)); ic = 1;
% fig = figure();
% % hold on
% c_map = parula(length(all_epsk)+2);
% 
% for epsilon_kappa = all_epsk
% 
%     % Setting for the material's time-modulation
%     Omega = 0.03; % modulation frequency
%     T = 2*pi/Omega;
%     phase_kappa = zeros(1,N); % modulation phases of kappa
%     phase_rho = zeros(1,N); % modulation phases of rho
%     for i = 1:(N-1)
%         phase_kappa(i+1) = pi/i;
%         phase_rho(i+1) = pi/i;
%     end
%     epsilon_rho = 0; % modulation amplitude of rho
%     rs = []; % Fourier coefficients of 1/rho
%     ks = []; % Fourier coefficients of 1/kappa
%     for j = 1:N
%         rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
%         ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
%         ks = [ks; ks_j];
%         rs = [rs; rs_j];
%     end
%     
%     % Calculate subwavelength resonant frequency
%     if N == 1
%         w_res = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0); % non-zero subwavelength resonant frequency
%     else
%         C = make_capacitance_finite(N,lij); % capacitance matrix
%         w_muller = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0); % subwavelength resonant frequencies
%         w_res = w_muller(real(w_muller)>=0); % positive subwavelength resonant frequencies
%     end
%     w_op = w_res(1)+0.0002; w0 = w_op;
%     k_op = w_op/v0; k0 = w0/v0;
%     
%     
%     % Define incident wave field  
%     vin_l = @(x,n) 1*exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % n-th mode of left incident wave
%     dx_vin_l = @(x,n) 1*sqrt(-1)*k0*exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % derivative of left incident wave
%     uin_r = @(x,t,n) 0; vin_r = @(x,n) 0; dx_vin_r = @(x,n) 0;
%     
%     % Build the linear system and solve it for the interior coefficients
%     MatcalA = getMatcalA(N,lij,xm,xp,k_tr,w_op,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
%     MatcalF = getF_lr(k_tr, N, delta, xm, xp, dx_vin_l, dx_vin_r); % vector \mathcal{F}
%     N_vin = getN_vin(k_tr, N, delta, xm, xp, w_op, Omega, v0, lij, vin_l, vin_r); % matrix-vector product \mathcal{N}v^{in}
%     RHS = MatcalF + N_vin;
%     sol = MatcalA\RHS; % solve for the interior coefficients, vector \mathbf{w}
%     
%     % Compute the transmission and reflection coefficients
%     [tns,rns] = get_tr(k_tr, w_op, Omega, rs, ks, vr, xm, xp, sol, lij, vin_l, vin_r);
%     O(:,ic) = abs(tns(:,1)).^2+abs(rns(:,end)).^2;
%     E(ic) = get_energy(tns(:,1),rns(:,end),w_op,Omega,v0,xm(1),xp(end),k_tr,t);
%     ic = ic+1;
% 
%     % Create plot
%     semilogy(-k_tr:k_tr,O(:,ic-1),'o','Color',c_map(ic-1,:),'DisplayName',strcat('$\varepsilon_{\kappa}=\,\,$',num2str(epsilon_kappa)),'MarkerSize',8,'LineWidth',2)
%     hold on
% 
% end
% 
% xlabel('$n$',Interpreter='latex',FontSize=18)
% ylabel('$|R_n|^2+|T_n|^2$',Interpreter='latex',FontSize=18)
% legend('show',interpreter='latex',fontsize=18)
% 
% figure()
% hold on
% plot(all_epsk,E,'o','MarkerSize',8,'LineWidth',2)
% xlabel('$\varepsilon_{\kappa}$',Interpreter='latex',FontSize=18)
% ylabel('Total Energy',Interpreter='latex',FontSize=18)


%% Test energy conservation over time

all_epsk = linspace(0,0.9,100);
% ts_cell = linspace(0,T,200);
% ts = [ts_cell,ts_cell+T];
ts = [0];
E = zeros(length(all_epsk),length(ts)); it = 1; ik = 1;
% fig1 = figure(1);
figure(1)
subplot(1,3,3)
hold on
c_map = parula(5);

for epsilon_kappa = all_epsk
    
    % Setting for the material's time-modulation
    Omega = 0.03; % modulation frequency
    T = 2*pi/Omega;
    phase_kappa = zeros(1,N); % modulation phases of kappa
    phase_rho = zeros(1,N); % modulation phases of rho
    for i = 1:(N-1)
        phase_kappa(i+1) = pi/i;
        phase_rho(i+1) = pi/i;
    end
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
    if N == 1
        w_res = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0); % non-zero subwavelength resonant frequency
    else
        C = make_capacitance_finite(N,lij); % capacitance matrix
        w_muller = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0); % subwavelength resonant frequencies
        w_res = w_muller(real(w_muller)>=0); % positive subwavelength resonant frequencies
    end
    w_op = w_res(1)+0.0002; w0 = w_op;
    k_op = w_op/v0; k0 = w0/v0;
    
    
    % Define incident wave field  
    vin_l = @(x,n) 1*exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % n-th mode of left incident wave
    dx_vin_l = @(x,n) 1*sqrt(-1)*k0*exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % derivative of left incident wave
    uin_r = @(x,t,n) 0; vin_r = @(x,n) 0; dx_vin_r = @(x,n) 0;
    
    % Build the linear system and solve it for the interior coefficients
    MatcalA = getMatcalA(N,lij,xm,xp,k_tr,w_op,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
    MatcalF = getF_lr(k_tr, N, delta, xm, xp, dx_vin_l, dx_vin_r); % vector \mathcal{F}
    N_vin = getN_vin(k_tr, N, delta, xm, xp, w_op, Omega, v0, lij, vin_l, vin_r); % matrix-vector product \mathcal{N}v^{in}
    RHS = MatcalF + N_vin;
    sol = MatcalA\RHS; % solve for the interior coefficients, vector \mathbf{w}

    it = 1;
    for t = ts
        % Compute the transmission and reflection coefficients
        [tns,rns] = get_tr(k_tr, w_op, Omega, rs, ks, vr, xm, xp, sol, lij, vin_l, vin_r);
        E(ik,it) = get_energy(tns(:,1),rns(:,end),w_op,Omega,v0,xm(1),xp(end),k_tr,t);
        it = it+1;
    end
    ik = ik + 1;

    % Plot results
%     plot(epsilon_kappa,E(ik-1,:),'*','Color',c_map(ik-1,:),'DisplayName',strcat('$\varepsilon_{\kappa}=\,\,$',num2str(epsilon_kappa)),'MarkerSize',8,'LineWidth',2)

end

% Plot results
plot(all_epsk,E,'-k','DisplayName',strcat('$N=\,\,$',num2str(N)),'MarkerSize',8,'LineWidth',2)
plot(linspace(0,0.9,50),ones(1,50),'--','Color',[.5 .5 .5],'HandleVisibility','off','MarkerSize',8,'LineWidth',2)
xlabel('$\varepsilon_{\kappa}$',Interpreter='latex',FontSize=18)
ylabel('Total Energy',Interpreter='latex',FontSize=18)
legend('show',interpreter='latex',fontsize=18)






