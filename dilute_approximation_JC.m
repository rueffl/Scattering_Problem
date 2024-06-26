% main script
format long
NN = 200;
xs = linspace(-4,4,NN);
uxs = zeros(NN,1);
uss = zeros(NN,1);

for i = 1:NN
    [uxs(i),uss(i)] = produce_plot(xs(i));
end
figure
plot(xs,real(uxs),'r');
figure
plot(xs,real(uss));

function [us,ux] = produce_plot(x)
    k_tr = 2; % truncation parameters as in remark 3.3
    N = 2; % number of the resonator
    li = ones(1,N)*0.005; % length of the resonators
    lij = ones(1,N-1); %lij(2:2:end) = lij(2:2:end)+1; % distance between the resonators
    L = sum(li)+sum(lij); % length of the unit cell
    Ls = zeros(2*N-1,1);
    Ls(1:2:end) = li;
    Ls(2:2:end) = lij;
    xipm = [0,cumsum(Ls)'];
    xm = xipm(1:2:end);
    xp = xipm(2:2:end);
    z = (xm+xp)./2;
    delta = 0.0001; % small contrast parameter
    vr = 1;
    v0 = 1;
    
    % Settings for modulation
    Omega = 0.03; % modulation frequency
    T = 2*pi/Omega;
    phase_kappa = zeros(1,N); % modulation phases of kappa
    phase_rho = zeros(1,N); % modulation phases of rho
    epsilon_kappa = 0; % modulation amplitude of kappa
    epsilon_rho = 0; % modulation amplitude of rho
    rs = [];
    ks = [];
    for j = 1:N
        rs_j = [epsilon_rho*exp(-1i*phase_rho(j))/2,1,epsilon_rho*exp(1i*phase_rho(j))/2];
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))/2,1,epsilon_kappa*exp(1i*phase_kappa(j))/2];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end
    
    % Find quasifrequencies
    w_muller = zeros(N,1);
    
    % Compute static case
    C = make_capacitance_finite(N,lij);
    w_static = get_capacitance_approx(0,0,li,Omega,phase_rho,phase_kappa,delta,C);
    w_out = get_capacitance_approx(epsilon_kappa,epsilon_rho,li,Omega,phase_rho,phase_kappa,delta,C);
    w_res = w_out(real(w_out)>=0);
    
    %green's function
    G = @(x,kn) exp(sqrt(-1)*kn*abs(x))./(2*sqrt(-1)*kn);
    
    t = 0;
    w_0 = mean(w_res);
    w = w_0 + 0.0001; %quasifrequency of incident wave
    k = w/v0;
    k_0 = w_0/v0; %wave number of incident wave
    A = getMatcalA(N, lij, xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0);
    F = getF(k_tr, N, delta, k, k_0, z);
    sol = linsolve(A,F);
    
    % Define coefficients alpha_n
    alphas = zeros(N,2*k_tr+1);
    all_cn = zeros(N,2*k_tr+1);
    gamma_s = zeros(N,2*k_tr+1);
    const = zeros(N,2*k_tr+1);
    for res = 1:N
        C = getC(k_tr, res, w, Omega, rs, ks, vr);
        [fi,lls] = eig(C,'vector');
        lambdas = zeros(1,2*k_tr+1);
        for j = -k_tr:k_tr
            lambdas(j+k_tr+1) = lls(k_tr-j+1);
        end
        for n = -k_tr:k_tr
            kn = (w+n*Omega)/v0;
            fn = zeros(1,2*k_tr+1);
            for j = -k_tr:k_tr
                lls = sqrt(lls);
                if k_tr-n+1 < 1 || k_tr-n+1 > 2*k_tr+1
                    fn(j+k_tr+1) = 0;
                else
                    fn(j+k_tr+1) = fi(k_tr-n+1,k_tr-j+1);
                end
            end
            
            alphas(res,n+k_tr+1) = get_cn(sol,kn,fn,lambdas,k_tr,z(res),res)*exp(-sqrt(-1)*k_0*z(res)); %We need alpha independent of res
            all_cn(res,n+k_tr+1) = get_cn(sol,kn,fn,lambdas,k_tr,z(res),res);
    
            gamma_s(res,n+k_tr+1) = operator_S(x, N, z, z, lij, k_tr, k, w, Omega, rs, ks, vr, sol, n)*(w-w_0);
            const(res,n+k_tr+1) = gamma_s(res,n+k_tr+1)./(G(x-z(res),(w+n*Omega)./v0));%.*exp(sqrt(-1).*k_0.*z(res))); %We need const=all_cn
    
        end
    end
    
    %incident wave
    % u_in = @(x,t) exp(sqrt(-1).*(k.*x+))
    
    ux = 0;
    for n = -k_tr:k_tr
        for j = 1:N
            ux = ux + all_cn(j,n+k_tr+1)*G(x-z(j),kn).*exp(sqrt(-1)*(w_0+n*Omega)*t);
        end
    end
    
    %compare with exact formula
    [us, xs] = prepare_plot(k_tr,li,lij,L,N,Omega,delta,phase_kappa,phase_rho,epsilon_kappa,epsilon_rho,xm,xp,vr,v0,[x]); %We need us = ux
end

