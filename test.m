clear
format long

% Settings for the material's structure
k_tr = 4; % truncation parameters as in remark 3.3
N = 2; % number of the resonator
spacing = 10; lij = ones(1,N-1).*spacing; % spacing between the resonators
len = 0.5; li = ones(1,N).*len; % length of the resonator
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
epsilon_kappa = 0.2; % modulation amplitude of kappa
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
    w_muller = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0,lij,xm,xp); % subwavelength resonant frequencies
    w_res = w_muller(real(w_muller)>=0); % positive subwavelength resonant frequencies
end
w_op = w_res(1)+ 0.0002; % operating frequency
all_w0 = linspace(0.001,4,100); % list of quasifrequencies of incident wave
k_op = w_op/v0; % operating wave number outside of the resonator
G = @(k,x) exp(sqrt(-1)*k*abs(x))./(2*sqrt(-1)*k); % Green's function
dx_G = @(k,x) x*exp(sqrt(-1)*k*abs(x))./(2*sqrt(-1)*abs(x)); % Derivative of Green's function

fig = figure();
hold on
c_map = parula(2*k_tr+1+2); % colors for plot

% iterate over w_0
for w0 = all_w0

    k0 = w0/v0; % wave number of incident frequency
    
    % Define incident wave field
    if N == 1
        vin_l = @(x,n) exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % n-th mode of left incident wave
        dx_vin_l = @(x,n) sqrt(-1)*k0*exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % derivative of left incident wave
        % uin_l = @(x,t,n) 0; vin_l = @(x,n) 0; dx_vin_l = @(x,n) 0;
        % uin_r = @(x,t,n) exp(sqrt(-1)*(-(k0).*x+w0.*t)).*(x>xp(end)).*(n==0); % right incident wave
        % vin_r = @(x,n) exp(sqrt(-1)*(-(k0).*x)).*(x>xp(end)).*(n==0); % n-th mode of right incident wave
        % dx_vin_r = @(x,n) -sqrt(-1)*k0*exp(sqrt(-1)*(-(k0).*x)).*(x>xp(end)).*(n==0); % derivative of right incident wave
        uin_r = @(x,t,n) 0; vin_r = @(x,n) 0; dx_vin_r = @(x,n) 0;

        % Build the linear system and solve it for the interior coefficients
        MatcalA = getMatcalA(N,lij,xm,xp,k_tr,w_op,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
        MatcalF = getF_lr(k_tr, N, delta, xm, xp, dx_vin_l, dx_vin_r); % vector \mathcal{F}
        sol = MatcalA\MatcalF; % solve for the interior coefficients, vector \mathbf{w}

        % Compute the frequency-scattering coefficients
        x = z(1) - 1;
        all_Ln = get_Lambdas_N1(x, k_tr, w_op, Omega, rs, ks, vr, v0, z, sol, vin_l, vin_r);

        % Create plot of the frequency-scattering coefficient
        for j = 1:2*k_tr+1
            if w0 == all_w0(1)
                plot(w0,all_Ln(j),'o','DisplayName',strcat('$\Lambda_',num2str(j-k_tr-1),'$'),'Color',c_map(j,:),'LineWidth',2,'MarkerSize',8)
            else
                plot(w0,all_Ln(j),'o','HandleVisibility','off','Color',c_map(j,:),'LineWidth',2,'MarkerSize',8)
            end
        end

    else

        all_Ln_l = zeros(N,2*k_tr+1); all_Ln_r = zeros(N,2*k_tr+1);
        for i = 1:N

            if i == 1
                vin_l = @(x,n) exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % n-th mode of left incident wave
                dx_vin_l = @(x,n) sqrt(-1)*k0*exp(sqrt(-1)*(k0).*x).*(x<xm(1)).*(n==0); % derivative of left incident wave
                % uin_l = @(x,t,n) 0; vin_l = @(x,n) 0; dx_vin_l = @(x,n) 0;
                % uin_r = @(x,t,n) exp(sqrt(-1)*(-(k0).*x+w0.*t)).*(x>xp(end)).*(n==0); % right incident wave
                % vin_r = @(x,n) exp(sqrt(-1)*(-(k0).*x)).*(x>xp(end)).*(n==0); % n-th mode of right incident wave
                % dx_vin_r = @(x,n) -sqrt(-1)*k0*exp(sqrt(-1)*(-(k0).*x)).*(x>xp(end)).*(n==0); % derivative of right incident wave
                uin_r = @(x,t,n) 0; vin_r = @(x,n) 0; dx_vin_r = @(x,n) 0;
            else
                vin_l = @(x,n) all_Ln_r(n+k_tr+1)*G((w_op+n*Omega)/vr,x-z(i))*(x<z(i)); % n-th mode of the left incident wave
                dx_vin_l = @(x,n) all_Ln_r(n+k_tr+1)*dx_G((w_op+n*Omega)/vr,x-z(i))*(x<z(i)); % derivative of left incident wave
                uin_r = @(x,t,n) 0; vin_r = @(x,n) 0; dx_vin_r = @(x,n) 0;
            end

            % Build the linear system and solve it for the interior coefficients
            MatcalA = getMatcalA(1,[],xm(i),xp(i),k_tr,w_op,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
            MatcalF = getF_lr(k_tr, 1, delta, xm(i), xp(i), dx_vin_l, dx_vin_r); % vector \mathcal{F}
            sol = MatcalA\MatcalF; % solve for the interior coefficients, vector \mathbf{w}

            % Compute the frequency-scattering coefficients
            xl = z(1) - 1;
            all_Ln_l(i,:) = get_Lambdas_N(xl, k_tr, w_op, Omega, rs, ks, vr, v0, z(i), sol, vin_l, vin_r);
            xr = z(1) + 1;
            all_Ln_r(i,:) = get_Lambdas_N(xr, k_tr, w_op, Omega, rs, ks, vr, v0, z(i), sol, vin_l, vin_r);

            % Create plot of the frequency-scattering coefficient
%             for j = 1:2*k_tr+1
%                 if w0 == all_w0(1)
%                     plot(w0,all_Ln_l(i,j),'o','DisplayName',strcat('$\Lambda_{',num2str(j-k_tr-1),',',num2str(i),'}$'),'Color',c_map(j,:),'LineWidth',2,'MarkerSize',8)
%                 else
%                     plot(w0,all_Ln_l(i,j),'o','HandleVisibility','off','Color',c_map(j,:),'LineWidth',2,'MarkerSize',8)
%                 end
%             end

        end

        % Create plot of the frequency-scattering coefficient
        for j = 1:2*k_tr+1
            if w0 == all_w0(1)
                plot(w0,all_Ln_l(:,j),'.','DisplayName',strcat('$\Lambda_{',num2str(j-k_tr-1),',i}$'),'Color',c_map(j,:),'LineWidth',2,'MarkerSize',8)
            else
                plot(w0,all_Ln_l(:,j),'.','HandleVisibility','off','Color',c_map(j,:),'LineWidth',2,'MarkerSize',8)
            end
        end

    end
    
    

end



xlabel('$\omega_{\mathrm{in}}$',Interpreter='latex',FontSize=18)
ylabel('$\Lambda_n$',Interpreter='latex',FontSize=18)
legend('show',interpreter='latex',fontsize=18)


