% main script
format long

%% Set parameters

% Settings for the structure
k_tr = 14; % truncation parameters as in remark 3.3
N = 6; % number of the resonator
li = ones(1,N); % length of the resonators
lij = ones(1,N-1); lij(2:2:end) = lij(2:2:end)+1; % distance between the resonators
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
for i = 1:(N-1)
    phase_kappa(i+1) = pi/i;
    phase_rho(i+1) = pi/i;
end
epsilon_kappa = 0.2; % modulation amplitude of kappa
epsilon_rho = 0; % modulation amplitude of rho
rs = [];
ks = [];
for j = 1:N
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j)),1,epsilon_rho*exp(-1i*phase_rho(j))];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j)),1,epsilon_kappa*exp(-1i*phase_kappa(j))];
    ks = [ks; ks_j];
    rs = [rs; rs_j];
end

% w_muller = zeros(1,N);
w_muller = zeros(1,N);
C = make_capacitance_finite(N,lij);
w_static = get_capacitance_approx(0,0,li,Omega,phase_rho,phase_kappa,delta,C);
% for i = 1:N
%     initial_guess = w_static(i);
%     w_muller(i) = muller(initial_guess,N,lij,xm,xp,k_tr,Omega,rs,ks,vr,delta,v0);
% %     [w_muller(i),p0] = search_vector_finite(initial_guess,N,lij,xm,xp,k_tr,Omega,rs,ks,vr,delta,v0);
% end
w_out = get_capacitance_approx_spec(epsilon_kappa,phase_kappa,Omega,delta,li,C);
w_muller = w_out(real(w_out)>=0);

w = w_muller(1) + 0.0001;
k = w/v0;
x = 6.5;

% gamma = zeros(2*k_tr+1,N);
gamma_s = zeros(2*k_tr+1,N);
for n = -k_tr:k_tr
    for j = 1:N

%         w0 = w_muller(j);
%         k_0 = w0/v0;
% 
%         A = getMatcalA(N, lij, xm, xp, k_tr, w0, Omega, rs, ks, vr, delta, v0);
%         F = getF(k_tr, N, delta, k, k_0, xm);
%         sol = linsolve(A,F);
% 
%         gamma(n+k_tr+1,j) = operator_S(x, N, xm, xp, lij, k_tr, k, w, Omega, rs, ks, vr, sol, n)*(w-w0);

        w0 = w_muller_s(j);
        k_0 = w0/v0;

        A = getMatcalA(N, lij, xm, xp, k_tr, w0, Omega, rs, ks, vr, delta, v0);
        F = getF(k_tr, N, delta, k, k_0, xm);
        sol = linsolve(A,F);

        gamma_s(n+k_tr+1,j) = operator_S(x, N, xm, xp, lij, k_tr, k, w, Omega, rs, ks, vr, sol, n)*(w-w0);

    end
end

gamma_approx = zeros(2*k_tr+1,N);
for n = -k_tr:k_tr
    for j = 1:N

        w0 = w_muller(j);
        zj = z(j);
        gamma_approx(n+k_tr+1,j) = exp(sqrt(-1)*(w0+n*Omega)*abs(x-zj));

    end
end

figure(1)
hold off
for j = 1:N

    subplot(2,3,j)
    plot(-k_tr:k_tr,real(gamma_s(:,j)),'r*')
    hold on
%     plot(-k_tr:k_tr,real(gamma(:,j)),'ro')
    plot(-k_tr:k_tr,real(gamma_approx(:,j)),'b*')
    if j == 1
        title([num2str(j), 'st resonator'], FontSize=40)
    elseif j == 2
        title([num2str(j), 'nd resonator'], FontSize=40)
    elseif j == 3
        title([num2str(j), 'rd resonator'], FontSize=40)
    else
        title([num2str(j), '-th resonator'], FontSize=40)
    end
    if j == N
        legend('Real Part', 'Imaginary Part', FontSize=36, Location='southoutside')
    end
    xlabel('$K$', Interpreter='latex', FontSize=36)

end




