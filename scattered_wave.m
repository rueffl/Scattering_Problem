
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

% Compute resonance frequencies
if N > 1
    C = make_capacitance_finite(N,lij);
    w_res = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0,lij,xm,xp);
    w_muller = w_res(real(w_res)>=0); 
else
    w_muller = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0);
    w_res = zeros(1,2); w_res(1,2) = w_muller; 
end

w = w_muller(1)+0.002; % operating frequency
win = 2.01*10^(-4); % incident frequency
k = w/v0; % mode 0 wave number outside of the resonators
k_in = win/v0; % wave number of the incident wave

% Define evaluation points
ls = 100;
us = zeros(ls,1);
xs = linspace(xm(1)-10, xp(N)+10,ls); % evaluation points

% Solve linear system for exterior coefficients
A = getMatcalA(N, lij, xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0);
F = getF(k_tr, N, delta, k, k_in, xm);
sol = A\F; % vector of solution coefficients inside the resonators

l = 1;
for x = xs % evaluate at each x
    for j = 1:N % iterate over the resonating frequencies
        wj = w_muller(j);
        for n = -k_tr:k_tr % iterate over the modes
            kn = (w+n*Omega)/v0; % mode n wave number outside the resonators
            knr = (w+n*Omega)/vr; % mode n wave number inside the resonators
            % check in which interval x lies
            if x <= xm(1) % x is before the first resonator, ie left of D_1
                % define beta_n^0
                ai = sol(1:2*N:end); % coefficients a_j^1
                bi = sol(2:2*N:end); % coefficients b_j^1
                lambda = ll(1, k_tr, w, Omega, rs, ks, vr); % sqrt of eigenvalues of matix C_1
                f = vs(1, k_tr, w, Omega, rs, ks, vr); % eigenvectors of matrix C_1
                f = f(n+k_tr+1,:);
                beta_n0 = 0;
                for c = -k_tr:k_tr
                    beta_n0 = beta_n0 + (ai(c+k_tr+1)*exp(sqrt(-1)*lambda(-c+k_tr+1)*xm(1))+bi(c+k_tr+1)*exp(-sqrt(-1)*lambda(-c+k_tr+1)*xm(1)))*f(-c+k_tr+1);
                end
                us(l) = us(l) + exp(sqrt(-1)*kn*xm(1))*beta_n0*exp(-sqrt(-1)*kn*x);
            elseif x >= xp(end) % x is after the last resonator, ie right of D_N
                % define alpha_n^N
                ai = sol(2*N-1:2*N:end); % coefficients a_j^N
                bi = sol(2*N:2*N:end); % coefficients b_j^N
                lambda = ll(N, k_tr, w, Omega, rs, ks, vr); % sqrt of eigenvalues of matix C_N
                f = vs(N, k_tr, w, Omega, rs, ks, vr); % eigenvectors of matrix C_N
                alpha_nN = 0;
                for c = -k_tr:k_tr
                    alpha_nN = alpha_nN + (ai(c+k_tr+1)*exp(sqrt(-1)*lambda(-c+k_tr+1)*xp(end))+bi(c+k_tr+1)*exp(-sqrt(-1)*lambda(-c+k_tr+1)*xp(end)))*f(-c+k_tr+1);
                end
                us(l) = us(l) + exp(-sqrt(-1)*kn*xp(end))*alpha_nN*exp(sqrt(-1)*kn*x);
            else
                for i = 1:N % check in which resonator or between which resonators x lies
                    if i < N
                        if xp(i) <= x && x <= xm(i+1) % between the i-th and (i+1)-th resonator
                            ab = operator_M(sol, i, k_tr, n, xm, xp, lij, k, w, Omega, rs, ks, vr);
                            us(l) = us(l) + (ab(1)*exp(sqrt(-1)*k*x)+ab(2)*exp(-sqrt(-1)*k*x))*exp(sqrt(-1)*(w_muller(j)+n*Omega)*t);
                        end
                        if xm(i) < x && x < xp(i) % inside the i-th resonator
                            lambda = ll(i, k_tr, w, Omega, rs, ks, vr); % sqrt of eigenvalues of matix C_i
                            f = vs(i, k_tr, w, Omega, rs, ks, vr); % eigenvectors of matrix C_i
                            ai = sol(2*i-1:2*N:end); % coefficients a_j^i
                            bi = sol(2*i:2*N:end); % coefficients b_j^i
                            s = 0;
                            for c = -k_tr:k_tr
                                s = s + (ai(c+k_tr+1)*exp(sqrt(-1)*lambda(-c+k_tr+1)*x)+bi(c+k_tr+1)*exp(-sqrt(-1)*lambda(-c+k_tr+1)*x))*f(-c+k_tr+1);
                            end
                            us(l) = us(l) + s;
                        end
                    end
                end
            end
        end
    end
    l = l+1;
end

% create plot
figure()
hold on
plot(xs,us,'-')





function lambdas = ll(i, k_tr, w, Omega, rs, ks, vr)
    C = getC(k_tr, i, w, Omega, rs, ks, vr);
    [~,lambdas] = eig(C,'vector');
    lambdas = sqrt(lambdas);
end

function fi = vs(i, k_tr, w, Omega, rs, ks, vr)
    C = getC(k_tr, i, w, Omega, rs, ks, vr);
    [fi,~] = eig(C,'vector');
end
