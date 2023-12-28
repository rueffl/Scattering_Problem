function [w_out] = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0)
% GET_CAPACITANCE_APPROX_HOT higher order capacitance approximation valid for any N>1 in 1D
%   epsilon_kappa:  modification amplitude of kappa
%   li:             length of the resonators
%   Omega:          modification frequency of kappa
%   phase_kappa:    phase of the modification of kappa
%   delta:          contrast parameter
%   C:              capacitance matrix
%   vr:             wave speed inside the resonators
%   v0:             wave speed outside the resonators

    T = 2*pi/Omega; % modulation period
    N = length(phase_kappa); % number of resonators
    invM2 = @(t) make_invM2(t,delta,vr,v0,li,epsilon_kappa,phase_kappa,Omega); % inverse of matrix M as a function of time
    M1 = @(t) make_M1(t,delta,vr,li,epsilon_kappa,phase_kappa,Omega,C); % matrix product of L and K'(t) multiplied by 1/(delta*vr^2)
    [w_out, cnd] = hill_exp(T,invM2,M1,N); % hill exponents of the ODE
    [w_out_real,order] = sort(real(w_out),'descend');
    w_out_imag = imag(w_out(order));
    w_out = w_out_real + sqrt(-1).*w_out_imag;

end

function out = make_invM2(t,delta,vr,v0,li,epsilon_kappa,phase_kappa,Omega)
% MAKE_INVM2 Creates the matrix M_2 and gives out its inverse
%   t:              time
%   delta:          contrast parameter
%   vr:             wave speed inside the resonators
%   v0:             wave speed outside the resonators
%   li:             length of the resonators
%   epsilon_kappa:  modulation amplitude of kappa
%   phase_kappa:    phase of the modification of kappa
%   Omega:          frequency of the modulation of kappa

    N = length(phase_kappa); % number of resonators
    invkappa = @(t) 1+epsilon_kappa.*cos(Omega.*t+phase_kappa); % function kappa(t)
    D = zeros(N,N); D(1,1) = 1; D(N,N) = 1; % matrix D
    L = diag(li); % matrix L
    Kt = diag(invkappa(t));
    Z = zeros(N,N);
    Id = eye(N);
    m2 = [Id, Z; D./v0, -L*Kt./(delta*vr^2)]; % matrix M_2
    out = inv(m2); % matrix M_2^{-1}

end

function out = make_M1(t,delta,vr,li,epsilon_kappa,phase_kappa,Omega,C)
% MAKE_LK Creates the matrix M_1
%   t:              time
%   delta:          contrast parameter
%   vr:             wave speed inside the resonators
%   v0:             wave speed outside the resonators
%   li:             length of the resonators
%   epsilon_kappa:  modulation amplitude of kappa
%   phase_kappa:    phase of the modification of kappa
%   Omega:          frequency of the modulation of kappa
%   C:              capacitance matrix

    N = length(phase_kappa); % number of resonators
    deriv_kappainv = @(t) -Omega.*epsilon_kappa.*sin(Omega.*t+phase_kappa); % function kappa'(t)
    LK = diag(li.*(1./(delta*vr^2)).*deriv_kappainv(t));
    Z = zeros(N,N);
    Id = eye(N);
    out = [Z, Id; C, LK];

end

function [w_out, cnd] = hill_exp(T,invM2,M1,N)
% HILL_EXP Computes the hill exponents of the 2-dimension 1st order ODE
%   T:      period of the time-modulation of kappa
%   invM2:  inverse of the matrix M_2, function of t
%   M1:     matrix M_1, function of t
%   N:      number of resonators

    W = zeros(2*N,2*N);
    II = eye(2*N,2*N);
    MM = @(t,y) (invM2(t)*M1(t))*y;
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