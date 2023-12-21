%% Set structure setting
clear 
format long

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

% Setting for the material's time-modulation
Omega = 0.03; % modulation frequency
T = 2*pi/Omega;
phase_kappa = zeros(1,N); % modulation phases of kappa
phase_rho = zeros(1,N); % modulation phases of rho
for i = 1:(N-1)
    phase_kappa(i+1) = pi;
    phase_rho(i+1) = pi;
end
epsilon_rho = 0; % modulation amplitude of rho


%% Apply vector search method

% tolerance for exiting the vector search method
tol = 1e-10;

% function to iterate with
p = @(ek) find_imagzero(ek,li,Omega,phase_kappa,delta,vr,v0,lij,T);

eks = linspace(0,3,800); p_eks = zeros(1,length(eks)); k = 1;
for ek = eks
    p_eks(k) = p(ek);
    k = k+1;
end
% figure();
% plot(eks,p_eks,'-','LineWidth',2);

% find minimum
[min_val,idx] = min(p_eks);
min_eks = eks(idx);

%% Compute energy for epsilon st Im(omega)=0

epsilon_kappa = min_eks; % modulation amplitude of kappa
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
w_op = w_res(1); w0 = w_op;
k_op = w_op/v0; k0 = k_op;

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

% Calculate energy
[tns,rns] = get_tr(k_tr, w_op, Omega, rs, ks, vr, xm, xp, sol, lij, vin_l, vin_r);
O = abs(tns(:,1)).^2+abs(rns(:,end)).^2;
% semilogy(-k_tr:k_tr,O,'*g','DisplayName',strcat('$\varepsilon_{\kappa}=\,\,$',num2str(epsilon_kappa)),'MarkerSize',8,'LineWidth',2)

%% Determine regions

w_res = @(epsilon_kappa) imag(omega_epsk(epsilon_kappa,li,Omega,phase_kappa,delta,vr,v0,lij,idx));
all_epsk = linspace(0,0.99,400);
w_epsk = zeros(1,length(all_epsk));

% prepare plot
% fig = figure();
% hold on
subplot(1,3,3)

for idx = 1:2*N

    change_idx = [];
    change_val = [];
    val_old = all_epsk(1);

    w_res = @(epsilon_kappa) imag(omega_epsk(epsilon_kappa,li,Omega,phase_kappa,delta,vr,v0,lij,idx));

    for i = 1:length(all_epsk)
        w_epsk(i) = w_res(all_epsk(i));
        % find change of sign
        if i > 1
            if w_epsk(i)*w_epsk(i-1) < 0
                change_idx = [change_idx,i];
            end
        end
        change_val = all_epsk(change_idx);
    end

    plot(all_epsk,w_epsk,'-k',markersize=5,linewidth=2)
    hold on
    xlabel('$\varepsilon_{\kappa}$','Interpreter','latex','fontsize',18)
    ylabel('$\mathrm{Im}(\omega)$','Interpreter','latex','fontsize',18)

end

%     % plot exact solutions
%     hold on
%     const = -2*1i*vr^2*delta/(li*v0);
%     w_epsk_exact = const./sqrt(1-all_epsk.^2);
%     plot(all_epsk,imag(w_epsk_exact),'r--',markersize=1,linewidth=2)



%% Functions

function [imag_w_res] = find_imagzero(epsilon_kappa,li,Omega,phase_kappa,delta,vr,v0,lij,T)
% FIND_IMAGZERO function whose root must be found in order to determine setting st imaginary part of resonant frequency is zero
%   epsilon_kappa:  time-modulation amplitude of kappa
%   li:             size of resonators
%   Omega:          modulation frequency of kappa
%   phase_kappa:    phase of time-modulated kappa
%   delta:          contrast parameter
%   vr:             wave speed inside resonators
%   v0:             wave speed in the background medium
%   lij:            spacing between the resonators

    N = length(phase_kappa);
    if N > 1
        C = make_capacitance_finite(N,lij); % capacitance matrix
        w_res = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0); % subwavelength resonant frequencies
%         w_res = w_res(real(w_res)>=0); % positive subwavelength resonant frequencies
    else
        w_res = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,li,delta,vr,v0); % non-zero subwavelength resonant frequency
%         w_res = get_omega_exact(delta,li,v0,vr,Omega,epsilon_kappa,phase_kappa,T);
    end
    imag_w_res = max(abs(imag(w_res))); % imaginary part of the swl resonant frequencies

end

function [w_res] = omega_epsk(epsilon_kappa,li,Omega,phase_kappa,delta,vr,v0,lij,idx)
% FIND_IMAGZERO function whose root must be found in order to determine setting st imaginary part of resonant frequency is zero
%   epsilon_kappa:  time-modulation amplitude of kappa
%   li:             size of resonators
%   Omega:          modulation frequency of kappa
%   phase_kappa:    phase of time-modulated kappa
%   delta:          contrast parameter
%   vr:             wave speed inside resonators
%   v0:             wave speed in the background medium
%   lij:            spacing between the resonators

    N = length(phase_kappa);
    if N > 1
        C = make_capacitance_finite(N,lij); % capacitance matrix
        w_res = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0); % subwavelength resonant frequencies
        [val,i] = sort(imag(w_res)); w_res = w_res(i); % sort omegas according to the real part
        [val,i] = max(imag(w_res));
    else
        w_res = [get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,li,delta,vr,v0),0]; % non-zero subwavelength resonant frequency
    end
    w_res = w_res(idx);

end





