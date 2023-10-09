function [w_out] = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0,lij,xm,xp)
% GET_CAPACITANCE_APPROX_HOT higher order capacitance approximation valid for any N>1 in 1D
%   epsilon_kappa:  modification amplitude of kappa
%   li:             length of the resonators
%   Omega:          modification frequency of kappa
%   phase_kappa:    phase of the modification of kappa
%   delta:          contrast parameter
%   C:              capacitance matrix
%   vr:             wave speed inside the resonators
%   v0:             wave speed outside the resonators
%   lij:            spacing between neighboring resonators
%   xm:             left boundary points
%   xp:             right boundary points

    T = 2*pi/Omega; % modulation period
    N = length(phase_kappa); % number of resonators
    invM = @(t) make_invM(t,delta,vr,v0,li,epsilon_kappa,phase_kappa,Omega,lij,xm,xp); % inverse of matrix M as a function of time
    LK = @(t) make_LK(t,delta,vr,li,epsilon_kappa,phase_kappa,Omega); % matrix product of L and K'(t) multiplied by 1/(delta*vr^2)
    [w_out, cnd] = hill_exp(T,invM,LK,N,C); % hill exponents of the ODE
    [w_out_real,order] = sort(real(w_out),'descend');
    w_out_imag = imag(w_out(order));
    w_out = w_out_real + sqrt(-1).*w_out_imag;

end

function out = make_invM(t,delta,vr,v0,li,epsilon_kappa,phase_kappa,Omega,lij,xm,xp)
% MAKE_INVM Creates the matrix M and gives out its inverse
%   t:              time
%   delta:          contrast parameter
%   vr:             wave speed inside the resonators
%   v0:             wave speed outside the resonators
%   li:             length of the resonators
%   epsilon_kappa:  modulation amplitude of kappa
%   phase_kappa:    phase of the modification of kappa
%   Omega:          frequency of the modulation of kappa
%   lij:            spacing between neighboring resonators
%   xm:             left boundary points
%   xp:             right boundary points

    kappa = @(t) 1./(1+epsilon_kappa.*cos(Omega.*t+phase_kappa)); % function kappa(t)
    T_offdiag_p = -xm(2:end).^2./lij(1:end);
    T_diag = [xm(1:end-1),0].*xp./[lij(1:end),1]+[0,xp(2:end)].*xm./[1,lij(1:end)];
    T_offdiag_m = -xp(1:end-1).^2./lij(1:end);
    T_mat = diag(T_offdiag_m,-1)+diag(T_diag)+diag(T_offdiag_p,1); % matrix T
    L = diag(li); % matrix L
    Kt = diag(kappa(t));
    m = (T_mat./(v0^2)+L*Kt./(delta*vr^2)); % matrix M
    out = inv(m); % matrix M^{-1}

end

function out = make_LK(t,delta,vr,li,epsilon_kappa,phase_kappa,Omega)
% MAKE_LK Creates the matrix product of L and K'(t) multiplied by 1/(vr^2*delta)
%   t:              time
%   delta:          contrast parameter
%   vr:             wave speed inside the resonators
%   v0:             wave speed outside the resonators
%   li:             length of the resonators
%   epsilon_kappa:  modulation amplitude of kappa
%   phase_kappa:    phase of the modification of kappa
%   Omega:          frequency of the modulation of kappa
%   lij:            spacing between neighboring resonators
%   xm:             left boundary points
%   xp:             right boundary points

    deriv_kappainv = @(t) -Omega.*epsilon_kappa.*sin(Omega.*t+phase_kappa); % function kappa'(t)
    out = diag(li.*(1./(delta*vr^2)).*deriv_kappainv(t));

end

function [w_out, cnd] = hill_exp(T,M,LK,N,C)
% HILL_EXP Computes the hill exponents of the 2-dimension 1st order ODE
%   T:      period of the time-modulation of kappa
%   M:      inverse of the matrix M
%   LK:     product of L an K'(t) multiplied by 1/(delta*vr^2)
%   N:      number of resonators
%   C:      capacitance matrix

    W = zeros(2*N,2*N);
    Z = zeros(N,N);
    I = eye(N,N);
    II = eye(2*N,2*N);
    MM = @(t,y) [Z, I; -M(t)*C, -M(t)*LK(t)]*y;
    for j = 1:2*N
        [~, wj] = ode45(MM,[0,T],II(:,j)); 
        W(:,j) = wj(end,:);
    end
    [U, D, V] = eig(W);
    w_out = (log(diag(D))/(1i*T));
    [out_real,ind] = sort(real(w_out),'descend');
    w_out = w_out(ind);
    cnd = cond(U);

end