%% Set structure setting
clear 
format long

% Settings for the material's structure
k_tr = 4; % truncation parameters as in remark 3.3
N = 6; % number of the resonator
% L = 2000; % length of the domain \mathcal{U}
% spacing = L/N; lij = ones(1,N-1).*spacing; % spacing between the resonators
spacing = 10; lij = ones(1,N-1).*spacing; % spacing between the resonators
len = 2; li = ones(1,N).*len; % length of the resonator
L = sum(li)+sum(lij); % length of the unit cell
Ls = zeros(2*N-1,1);
Ls(1:2:end) = li;
Ls(2:2:end) = lij;
xipm = [0,cumsum(Ls)']-(len*(N/2)+spacing*(N-1)/2); % all boundary points, make sure the resonators are aligned symmetrically wrt 0
xm = xipm(1:2:end); % LHS boundary points
xp = xipm(2:2:end); % RHS boundary points
z = (xm+xp)./2; % centers of resonators
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
if N > 1
    w_res = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0)*ones(1,N);
    w_op = w_res(1)+ 0.0002; % operating frequency
    w0 = w_op; %w0 = real(w_res(1))+0.02; % quasifrequency of incident wave
else
    w_res = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0); % non-zero subwavelength resonant frequency
    w_op = w_res+0.0002; % operating frequency
    w0 = w_op; % quasifrequency of incident wave
end
k_op = w_op/v0; % operating wave number outside of the resonator
k0 = w0/v0; % wave number of incident frequency

% Define relevant functions
uin = @(x,t) 1*exp(sqrt(-1)*((k0).*x+w0.*t)).*(x<xm(1)); % incident wave
vin = @(x,n) 1*exp(sqrt(-1)*k0*x)*(n==0)*(x<xm(1));
dx_vin = @(x,n) 1*sqrt(-1)*k0*exp(sqrt(-1)*k0*x)*(n==0)*(x<xm(1));
G = @(k,x) exp(sqrt(-1)*k*abs(x))./(2*sqrt(-1)*k); % Green's function

% Define evaluation points
len_xs = 800;
len_zs = 80;
xs = linspace(xm(1)-spacing,xp(end)+spacing,len_xs);
zs = zeros(N,len_zs);
for i = 1:N
    zs(i,:) = linspace(xm(i),xp(i),len_zs);
end


%% Compute the scattered wave using formula (31)

% Compute solution coefficients
MatcalA = getMatcalA(N,lij,xm,xp,k_tr,w_op,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
MatcalF = getF_lr(k_tr, N, delta, xm, xp, dx_vin, @(x,n) 0); % vector \mathcal{F}
N_vin = getN_vin(k_tr, N, delta, xm, xp, w_op, Omega, v0, lij, vin, @(x,n) 0); % matrix-vector product \mathcal{N}v^{in}
RHS = MatcalF + N_vin;
sol = MatcalA\RHS; % solve for the interior coefficients, vector \mathbf{w}

us_eval_x = zeros(1,len_xs);
us_eval_z = zeros(1,len_zs);
usx = @(x) uin(x,t) + get_us(x, t, N, xm, xp, lij, k_tr, v0, w_op, Omega, rs, ks, vr, sol, w_res, k0, vin); % scattered wave field as a function of x for fixed time t, according to formula (31)

for i = 1:len_xs
    us_eval_x(i) = usx(xs(i));
end

% Create plot
fig = figure();
fig.Position = [263,725,982,352];
subplot(1,2,1)
set(gca,'FontSize',14)
hold on
plot(xs,real(us_eval_x),'-','DisplayName','Exact',markersize=8,linewidth=2)
xlim([xs(1) xs(end)])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel(strcat('Re$(u^{\mathrm{sc}}(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)

subplot(1,2,2)
set(gca,'FontSize',14)
hold on
plot(xs,imag(us_eval_x),'-','DisplayName','Exact',markersize=8,linewidth=2)
xlim([xs(1) xs(end)])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel(strcat('Im$(u^{\mathrm{sc}}(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)

for i = 1:N
    for j = 1:len_zs
        us_eval_z(j) = usx(zs(i,j));
    end
    subplot(1,2,1)
    plot(zs(i,:),real(us_eval_z),'r-','DisplayName','Exact',markersize=8,linewidth=2)
    subplot(1,2,2)
    plot(zs(i,:),imag(us_eval_z),'r-','DisplayName','Exact',markersize=8,linewidth=2)
end


%% Iterate over epsilon and create plot

all_epsk = [0,0.1,0.2,0.3,0.4];

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
z = (xm+xp)./2; % centers of resonators
delta = 0.0001; % small contrast parameter
t = 0; % time

vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators

% Define evaluation points
len_xs = 800;
len_zs = 80;
xs = linspace(xm(1),xp(end)+spacing,len_xs);
zs = zeros(N,len_zs);
for i = 1:N
    zs(i,:) = linspace(xm(i),xp(i),len_zs);
end

% Setting for the material's time-modulation
Omega = 0.03; % modulation frequency
T = 2*pi/Omega;
phase_kappa = zeros(1,N); % modulation phases of kappa
phase_rho = zeros(1,N); % modulation phases of rho
for i = 1:(N-1)
    phase_kappa(i+1) = pi/i;
    phase_rho(i+1) = pi/i;
end

% prepare figure
fig = figure();
fig.Position = [263,725,982,352];
c_map = parula(length(all_epsk)+2); ic = 1;
O = zeros(2*k_tr+1,length(all_epsk));

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
    C = make_capacitance_finite(N,lij); % capacitance matrix
    w_muller = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0); % subwavelength resonant frequencies
    w_res = w_muller(real(w_muller)>=0); % positive subwavelength resonant frequencies
    w_op = w_res(1)+ 0.0002; % operating frequency
    w0 = w_op;%real(w_res(1))+0.02; % quasifrequency of incident wave
    k_op = w_op/v0; % operating wave number outside of the resonator
    k0 = w0/v0; % wave number of incident frequency

    % Compute solution coefficients
    MatcalA = getMatcalA(N,lij,xm,xp,k_tr,w_op,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
    MatcalF = getF_lr(k_tr, N, delta, xm, xp, dx_vin, @(x,n) 0); % vector \mathcal{F}
    N_vin = getN_vin(k_tr, N, delta, xm, xp, w_op, Omega, v0, lij, vin, @(x,n) 0); % matrix-vector product \mathcal{N}v^{in}
    RHS = MatcalF + N_vin;
    sol = MatcalA\RHS; % solve for the interior coefficients, vector \mathbf{w}
    
    us_eval_x = zeros(1,len_xs);
    us_eval_z = zeros(1,len_zs);
    usx = @(x) uin(x,t) + get_us(x, t, N, xm, xp, lij, k_tr, v0, w_op, Omega, rs, ks, vr, sol, w_res, k0, vin); % scattered wave field as a function of x for fixed time t, according to formula (31)

    for i = 1:len_xs
        us_eval_x(i) = usx(xs(i));
    end
    
    % plot results
    subplot(1,2,1)
    set(gca,'FontSize',14)
    hold on
    plot(xs,real(us_eval_x),'-','DisplayName',strcat('$\varepsilon_{\kappa}=$',num2str(epsilon_kappa)),'Color',c_map(ic,:),markersize=8,linewidth=2)
    
    subplot(1,2,2)
    set(gca,'FontSize',14)
    hold on
    plot(xs,imag(us_eval_x),'-','DisplayName',strcat('$\varepsilon_{\kappa}=$',num2str(epsilon_kappa)),'Color',c_map(ic,:),markersize=8,linewidth=2)
    
    for i = 1:N
        for j = 1:len_zs
            us_eval_z(j) = usx(zs(i,j));
        end
        subplot(1,2,1)
        plot(zs(i,:),real(us_eval_z),'-','HandleVisibility','off','Color',c_map(ic,:),markersize=8,linewidth=2)
        subplot(1,2,2)
        plot(zs(i,:),imag(us_eval_z),'-','HandleVisibility','off','Color',c_map(ic,:),markersize=8,linewidth=2)
    end

    [tns,rns] = get_tr(k_tr, w_op, Omega, rs, ks, vr, xm, xp, sol, lij, vin, @(x,n) 0);
    O(:,ic) = abs(tns(:,1)).^2+abs(rns(:,end)).^2;
    ic = ic+1;
    
end

subplot(1,2,1)
xlim([xs(1) xs(end)])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel(strcat('Re$(u(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)

subplot(1,2,2)
xlim([xs(1) xs(end)])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel(strcat('Im$(u(x,$',num2str(t),'$))$'),Interpreter='latex',FontSize=18)

legend('show',interpreter='latex',fontsize=18,location='southoutside')


%% Iterate over time and create plot

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
z = (xm+xp)./2; % centers of resonators
delta = 0.0001; % small contrast parameter
ts = linspace(0,100,11); % time

vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators

% Define evaluation points
len_xs = 800;
len_zs = 80;
xs = linspace(xm(1),xp(end)+spacing,len_xs);
zs = zeros(N,len_zs);
for i = 1:N
    zs(i,:) = linspace(xm(i),xp(i),len_zs);
end

% Setting for the material's time-modulation
Omega = 0.03; % modulation frequency
T = 2*pi/Omega;
phase_kappa = zeros(1,N); % modulation phases of kappa
phase_rho = zeros(1,N); % modulation phases of rho
for i = 1:(N-1)
    phase_kappa(i+1) = pi/i;
    phase_rho(i+1) = pi/i;
end
epsilon_kappa = 0.4;
epsilon_rho = 0;
rs = []; % Fourier coefficients of 1/rho
ks = []; % Fourier coefficients of 1/kappa
for j = 1:N
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
    ks = [ks; ks_j];
    rs = [rs; rs_j];
end

% prepare figure
fig = figure();
fig.Position = [2226, 623, 1219, 345];
c_map = parula(length(ts)+2); ic = 1;

for t = ts
    
    % Calculate subwavelength resonant frequency
    C = make_capacitance_finite(N,lij); % capacitance matrix
    w_muller = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0); % subwavelength resonant frequencies
    w_res = w_muller(real(w_muller)>=0); % positive subwavelength resonant frequencies
    w_op = w_res(1)+ 0.0002; % operating frequency
    w0 = w_op;%real(w_res(1))+0.02; % quasifrequency of incident wave
    k_op = w_op/v0; % operating wave number outside of the resonator
    k0 = w0/v0; % wave number of incident frequency

    % Define relevant functions
    uin = @(x,t) 1*exp(sqrt(-1)*((k0).*x+w0.*t)).*(x<xm(1)); % incident wave
    vin = @(x,n) 1*exp(sqrt(-1)*k0*x)*(n==0)*(x<xm(1));
    dx_vin = @(x,n) 1*sqrt(-1)*k0*exp(sqrt(-1)*k0*x)*(n==0)*(x<xm(1));

    % Compute solution coefficients
    MatcalA = getMatcalA(N,lij,xm,xp,k_tr,w_op,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
    MatcalF = getF_lr(k_tr, N, delta, xm, xp, dx_vin, @(x,n) 0); % vector \mathcal{F}
    N_vin = getN_vin(k_tr, N, delta, xm, xp, w_op, Omega, v0, lij, vin, @(x,n) 0); % matrix-vector product \mathcal{N}v^{in}
    RHS = MatcalF + N_vin;
    sol = MatcalA\RHS; % solve for the interior coefficients, vector \mathbf{w}
    
    us_eval_x = zeros(1,len_xs);
    us_eval_z = zeros(1,len_zs);
    usx = @(x) uin(x,t) + get_us(x, t, N, xm, xp, lij, k_tr, v0, w_op, Omega, rs, ks, vr, sol, w_res, k0, vin); % scattered wave field as a function of x for fixed time t, according to formula (31)

    for i = 1:len_xs
        us_eval_x(i) = usx(xs(i));
    end
    
    % plot results
    subplot(1,2,1)
    set(gca,'FontSize',14)
    hold on
    plot(xs,real(us_eval_x),'-','DisplayName',strcat('$t=$',num2str(t)),'Color',c_map(ic,:),markersize=8,linewidth=2)
    
    subplot(1,2,2)
    set(gca,'FontSize',14)
    hold on
    plot(xs,imag(us_eval_x),'-','DisplayName',strcat('$t=$',num2str(t)),'Color',c_map(ic,:),markersize=8,linewidth=2)
    
    for i = 1:N
        for j = 1:len_zs
            us_eval_z(j) = usx(zs(i,j));
        end
        subplot(1,2,1)
        plot(zs(i,:),real(us_eval_z),'-','HandleVisibility','off','Color',c_map(ic,:),markersize=8,linewidth=2)
        subplot(1,2,2)
        plot(zs(i,:),imag(us_eval_z),'-','HandleVisibility','off','Color',c_map(ic,:),markersize=8,linewidth=2)
    end

    ic = ic+1;
    
end

subplot(1,2,1)
xlim([xs(1) xs(end)])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel('Re$(u(x,t))$',Interpreter='latex',FontSize=18)

subplot(1,2,2)
xlim([xs(1) xs(end)])
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel('Im$(u(x,t))$',Interpreter='latex',FontSize=18)

legend('show',interpreter='latex',fontsize=18,location='southoutside')

%% Create 3D plot of the total wave field

% Settings for the material's structure
k_tr = 4; % truncation parameters as in remark 3.3
N = 6; % number of the resonator
spacing = 10; lij = ones(1,N-1).*spacing; % spacing between the resonators
len = 2; li = ones(1,N).*len; % length of the resonator
Ls = zeros(2*N-1,1);
Ls(1:2:end) = li;
Ls(2:2:end) = lij;
xipm = [0,cumsum(Ls)']-(len*(N/2)+spacing*(N-1)/2); % all boundary points, make sure the resonators are aligned symmetrically wrt 0
xm = xipm(1:2:end); % LHS boundary points
xp = xipm(2:2:end); % RHS boundary points

xs = linspace(xm(1),xp(end),200);
ts = linspace(0,220,200);
uxt = zeros(200,200);

for it = 1:200
    for ix = 1:200
        uxt(ix,it) = total_u_xt(xs(ix),ts(it));
    end
end

s = surf(xs,ts,abs(uxt),'EdgeColor','interp');
s.EdgeColor = 'none';
xlabel('$x$',interpreter='latex')
ylabel('$t$',interpreter='latex')
zlabel('$u(x,t)$',interpreter='latex')


%% function for 3D plot

function uxt = total_u_xt(xi,ti)

    % Settings for the material's structure
    k_tr = 4; % truncation parameters as in remark 3.3
    N = 6; % number of the resonator
    spacing = 10; lij = ones(1,N-1).*spacing; % spacing between the resonators
    len = 2; li = ones(1,N).*len; % length of the resonator
    Ls = zeros(2*N-1,1);
    Ls(1:2:end) = li;
    Ls(2:2:end) = lij;
    xipm = [0,cumsum(Ls)']-(len*(N/2)+spacing*(N-1)/2); % all boundary points, make sure the resonators are aligned symmetrically wrt 0
    xm = xipm(1:2:end); % LHS boundary points
    xp = xipm(2:2:end); % RHS boundary points
    delta = 0.0001; % small contrast parameter
    
    vr = 1; % wave speed inside the resonators
    v0 = 1; % wave speed outside the resonators
    
    % Setting for the material's time-modulation
    Omega = 0.03; % modulation frequency
    phase_kappa = zeros(1,N); % modulation phases of kappa
    phase_rho = zeros(1,N); % modulation phases of rho
    for i = 1:(N-1)
        phase_kappa(i+1) = pi/i;
        phase_rho(i+1) = pi/i;
    end
    epsilon_kappa = 0.2;
    epsilon_rho = 0;
    rs = []; % Fourier coefficients of 1/rho
    ks = []; % Fourier coefficients of 1/kappa
    for j = 1:N
        rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end
        
    % Calculate subwavelength resonant frequency
    C = make_capacitance_finite(N,lij); % capacitance matrix
    w_muller = get_capacitance_approx_hot(epsilon_kappa, li, Omega, phase_kappa, delta, C, vr, v0); % subwavelength resonant frequencies
    w_res = w_muller(real(w_muller)>=0); % positive subwavelength resonant frequencies
    w_op = w_res(1)+ 0.0002; % operating frequency
    w0 = w_op;%real(w_res(1))+0.02; % quasifrequency of incident wave
    k0 = w0/v0; % wave number of incident frequency

    % Define relevant functions
    uin = @(x,t) 1*exp(sqrt(-1)*((k0).*x+w0.*t)).*(x<xm(1)); % incident wave
    vin = @(x,n) 1*exp(sqrt(-1)*k0*x)*(n==0)*(x<xm(1));
    dx_vin = @(x,n) 1*sqrt(-1)*k0*exp(sqrt(-1)*k0*x)*(n==0)*(x<xm(1));

    % Compute solution coefficients
    MatcalA = getMatcalA(N, lij, xm, xp, k_tr, w_op, Omega, rs, ks, vr, delta, v0); % matrix \mathcal{A}
    MatcalF = getF_lr(k_tr, N, delta, xm, xp, dx_vin, @(x,n) 0); % vector \mathcal{F}
    N_vin = getN_vin(k_tr, N, delta, xm, xp, w_op, Omega, v0, lij, vin, @(x,n) 0); % matrix-vector product \mathcal{N}v^{in}
    RHS = MatcalF + N_vin;
    sol = MatcalA\RHS; % solve for the interior coefficients, vector \mathbf{w}

    uxt = uin(xi,ti) + get_us(xi, ti, N, xm, xp, lij, k_tr, v0, w_op, Omega, rs, ks, vr, sol, w_res, k0, vin); % scattered wave field as a function of x for fixed time t, according to formula (31)

end


