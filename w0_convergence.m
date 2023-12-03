format long

% Settings for the structure
k_tr = 4; % truncation parameters as in remark 3.3
N = 6; % number of the resonator
li = ones(1,N).*0.025; % length of the resonators
lij_s = linspace(500,500000,100); % distance between the resonators
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


all_w = zeros(length(lij_s),2*N);
j = 1;
fig = figure();
hold on
gray = [.7 .7 .7];
for len = lij_s
    lij = ones(1,N-1).*len;
    L = sum(li)+sum(lij); % length of the unit cell
    Ls = zeros(2*N-1,1);
    Ls(1:2:end) = li;
    Ls(2:2:end) = lij;
    xipm = [0,cumsum(Ls)']; % all boundary points
    xm = xipm(1:2:end); % LHS boundary points
    xp = xipm(2:2:end); % RHS boundary points
    C = make_capacitance_finite(N,lij);
    w_out = get_capacitance_approx(epsilon_kappa,epsilon_rho,li,Omega,phase_rho,phase_kappa,delta,C);
    w_out_hot = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0); % subwavelength resonant frequencies; get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0,lij,xm,xp);
    if len == lij_s(end)
        plot(ones(1,2*N).*len,real(w_out),'.','Color',gray,'MarkerSize',8,'DisplayName','Order $O(\delta)$ approximation')
        plot(ones(1,2*N).*len,real(w_out_hot),'k.','MarkerSize',8,'DisplayName','Order $O(\delta^{2/3})$ approximation')
    else
        plot(ones(1,2*N).*len,real(w_out),'.','Color',gray,'MarkerSize',8,'HandleVisibility','off')
        plot(ones(1,2*N).*len,real(w_out_hot),'k.','MarkerSize',8,'HandleVisibility','off')
    end
    all_w(j,:) = w_out;
    j = j+1;
end

% w1 = -sqrt(-1)*vr*log(1+2*vr*delta/(v0-vr*delta));
w1 = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,lij_s(1),delta,vr,v0);
plot(lij_s(end),real(w1),'r*','DisplayName','$\omega_0$ for $N=1$')

fig.Position = [2277,239,1699,685];
xlabel('$\ell_{i(i+1)}$',interpreter='latex',fontsize=20)
ylabel('$\omega_j$',interpreter='latex',fontsize=20)
legend('show',interpreter='latex',fontsize=18,location='best')

%% N = 1

format long

% Settings for the structure
k_tr = 2; % truncation parameters as in remark 3.3
N = 1; % number of the resonator
lij = ones(1,N-1); % spacing between the resonators
len = 1;
li = ones(1,N).*len; % length of the resonator
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

guess = zeros(1,2);
guess(2) = -(v0*vr*log((delta^2*vr^2 + 2*delta*v0*vr + v0^2)/(v0 - delta*vr)^2)*sqrt(-1))/(li(1)*(v0 + vr));%-sqrt(-1)*vr*log(1+2*vr*delta/(v0-vr*delta));
w_out = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0);
w_res = zeros(1,2);
for k = 1:2*N
    w_res(k) = muller(guess(k),N,lij,xm,xp,k_tr,Omega,rs,ks,vr,delta,v0);
end
approx_w(j,2) = -(v0*vr*log((delta^2*vr^2 + 2*delta*v0*vr + v0^2)/(v0 - delta*vr)^2)*sqrt(-1))/(li(1)*(v0 + vr));





%% N=1, static case, iterate over li

format long

% Settings for the structure
k_tr = 2; % truncation parameters as in remark 3.3
N = 1; % number of the resonator
lij = ones(1,N-1); % spacing between the resonators
li_s = exp(linspace(log(0.05),log(10),40)); % distance between the resonators
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


w_res = zeros(length(li_s),2*N);
w_approx = zeros(length(li_s),2*N);
j = 1;

for len = li_s
    li = ones(1,N).*len;
    xm = 0; xp = xm+li(1);
    w_out = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0);
    w_res(j,2) = w_out;
    for k = 1:2*N
        w_approx(j,k) = muller(w_res(j,k)+0.01,N,lij,xm,xp,k_tr,Omega,rs,ks,vr,delta,v0);
    end
    j = j+1;
end

figure()
hold on
plot(li_s,imag(w_res(:,1)),'k*','DisplayName','Approximation',markersize=8)
plot(li_s,imag(w_res(:,2)),'k*',markersize=8)
plot(li_s,imag(w_approx(:,1)),'r.','DisplayName','Muller`s Method',markersize=8)
plot(li_s,imag(w_approx(:,2)),'r.',markersize=8)
legend('show',interpreter='latex',fontsize=18,location='southoutside')
% legend('Approximation','Im$(\omega_i(\delta))$',interpreter='latex',fontsize=24)
xlabel('$\ell_{i}$',interpreter='latex',fontsize=20)
ylabel('Im$(\omega_j)$',interpreter='latex',fontsize=20)

figure()
hold on
plot(li_s,real(w_res(:,1)),'k*','DisplayName','Approximation',markersize=8)
plot(li_s,real(w_res(:,2)),'k*',markersize=8)
plot(li_s,real(w_approx(:,1)),'r.','DisplayName','Muller`s Method',markersize=8)
plot(li_s,real(w_approx(:,2)),'r.',markersize=8)
legend('show',interpreter='latex',fontsize=18,location='southoutside')
% legend('Approximation','Re$(\omega_i(\delta))$',interpreter='latex',fontsize=24)
xlabel('$\ell_{i}$',interpreter='latex',fontsize=20)
ylabel('Re$(\omega_j)$',interpreter='latex',fontsize=20)

% compute error between approximation and Muller's method
abs_error = abs(w_res-w_approx);


%% N=1, iterate over epsilon_kappa

format long

% Settings for the structure
k_tr = 2; % truncation parameters as in remark 3.3
N = 1; % number of the resonator
lij = ones(1,N-1); % spacing between the resonators
% li_s = linspace(0.05,20,80); % distance between the resonators
li_s = exp(linspace(log(0.05),log(10),40));
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
eks = linspace(0,0.8,9); % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho

figure()
c_map = parula(length(eks)+2); iR = 1;
for epsilon_kappa = eks 

    rs = [];
    ks = [];
    for j = 1:N
        rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end
    
    w_res = zeros(length(li_s),2*N);
    j = 1;
    
    for len = li_s
        li = ones(1,N).*len;
        xm = 0; xp = xm+li(1);
        w_out = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0);
        w_res(j,2) = w_out;
        j = j+1;
    end
    
    loglog(li_s,imag(w_res(:,2)),'-','Color',c_map(iR,:),markersize=8,linewidth=2)
    hold on
    iR = iR+1;

end

legend('$\epsilon_{\kappa}=0$','$\epsilon_{\kappa}=0.1$','$\epsilon_{\kappa}=0.2$',...
    '$\epsilon_{\kappa}=0.3$','$\epsilon_{\kappa}=0.4$','$\epsilon_{\kappa}=0.5$',...
    '$\epsilon_{\kappa}=0.6$','$\epsilon_{\kappa}=0.7$','$\epsilon_{\kappa}=0.8$',...
    '$\epsilon_{\kappa}=0.9$','$\epsilon_{\kappa}=1.0$','Approximation',...
    interpreter='latex',fontsize=24,location='east outside')

xlabel('$\ell_{i}$',interpreter='latex',fontsize=20)
ylabel('Im$(\omega_j)$',interpreter='latex',fontsize=20)




