%% Check if formula (53) agrees with exact computations in the dilute regime

% Settings for the structure
k_tr = 2; % truncation parameters as in remark 3.3
N = 6; % number of the resonator
spacing = 1; lij = ones(1,N-1).*spacing; % spacing between the resonators
len = 0.1; li = ones(1,N).*len; % length of the resonator
L = sum(li)+sum(lij); % length of the unit cell
Ls = zeros(2*N-1,1);
Ls(1:2:end) = li;
Ls(2:2:end) = lij;
xipm = [0,cumsum(Ls)']; % all boundary points
xm = xipm(1:2:end); % LHS boundary points
xp = xipm(2:2:end); % RHS boundary points
z = (xm+xp)./2; % centers of resonators
delta = 0.0001; % small contrast parameter
t = 0; % time

vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators

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
rs = []; % Fourier coefficients of 1/rho
ks = []; % Fourier coefficients of 1/kappa
for j = 1:N
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
    ks = [ks; ks_j];
    rs = [rs; rs_j];
end

% exact computation
[us, xs] = prepare_plot(k_tr,li,lij,L,N,Omega,delta,phase_kappa,phase_rho,epsilon_kappa,epsilon_rho,xm,xp,vr,v0,[]);

% use formula (53) 
% consider the single-resonator setting to compute \Lambda_n for all n=-K,...,K
w_out = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0); % compute \omega_1
w_res = zeros(1,2); w_res(1,2) = w_out; 
w0 = w_out + 0.0002; % operating frequency
w = w0 + 0.00000001; % quasifrequency of incident wave
k = w/v0; % wave number of incident wave
k_0 = w0/v0; % wave number of operating wave

all_Ln = get_Lambdas(xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0, k, k_0)'; % \Lambda_n for all n=-K,...,K
G = @(k,x) exp(sqrt(-1).*k.*abs(x))./(2.*sqrt(-1).*k); % Green's function 
uin = @(x,t) exp(k.*x+w.*t); % incident wave

ls = 100;
sumn = zeros(ls,1);
xs = linspace(xm(1)-1, xp(N)+1,ls); % evaluation points

for l = 1:ls
    x = xs(l); % evaluation point
    for n = -k_tr:k_tr % iterate over modes
        temp_sum = 0;
        kn = (w0+n*Omega)./v0;
        for j = 1:N % iterate over resonators
            temp_sum = temp_sum + G(kn,x-z(j)).*uin(z(j),t);
        end
        sumn(l) = sumn(l) + all_Ln(n+k_tr+1)*temp_sum*exp(sqrt(-1)*n*Omega*t);
    end
    sumn(l) = sumn(l)./(w0-w_out);
end

% create plot
figure()
hold on
plot(xs,real(sumn),'-','DisplayName','Formula (52)',markersize=8,linewidth=2)
plot(xs,real(us),'-','DisplayName','Exact',markersize=8,linewidth=2)
legend('show',Interpreter='latex',FontSize=18,location='southoutside')
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel('Real Part',Interpreter='latex',FontSize=18)

figure()
hold on
plot(xs,imag(sumn),'-','DisplayName','Formula (52)',markersize=8,linewidth=2)
plot(xs,imag(us),'-','DisplayName','Exact',markersize=8,linewidth=2)
legend('show',Interpreter='latex',FontSize=18,location='southoutside')
xlabel('$x$',Interpreter='latex',FontSize=18)
ylabel('Imaginary Part',Interpreter='latex',FontSize=18)



