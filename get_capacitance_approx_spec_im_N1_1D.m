function w_out = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,phase_kappa,Omega,l,delta,vr,v0)
%GET_CAPACITANCE_APPROX_SPEC_IM_N1_1D Computes the subwavelength
%quasifrequencies of a single-resonator system in 1D
%   epsilon_kappa:  modulation amplitude of kappa
%   phase_kappa:    phase shifts of kappa
%   Omega:          frequency of kappa
%   l:              length of the resonator
%   delta:          contrast parameter
%   vr:             wave speed inside the resonator
%   v0:             wave speed outside the resonator

    M = 1; % Number of Fourier coefficients of 1/\kappa
    N = 1;
    N_fourier = 4; % Length of Fourier series approx
    
    K_mod = zeros(2*M+1,N);
    for i = 1:N
        K_mod(:,i) = [epsilon_kappa/2*exp(-1i*phase_kappa(i)); 1; epsilon_kappa/2*exp(1i*phase_kappa(i))]; % Fourier coefficients of 1/\kappa
    end

    ns = -N_fourier:N_fourier;

    NN = 2*N_fourier+1;
    O = diag(ns.'*Omega);
    e = ones(NN,1);
    K = zeros(NN);
    for m = -M:M
        K = K+diag(e(1:NN-abs(m))*K_mod(m+M+1,i),m);
    end
    iK = inv(K); %% Fourier coefficients of \kappa

    c = 2*delta*(vr)^2/(v0*l);
    mat = -O-1i*c.*iK; 

    w_out = sort(eigs(mat,2*N,'smallestabs'),'ComparisonMethod','real'); % The eigenvalues of "mat" are approximately \omega + n\Omega for |n| < N_fouier. Taking the smallest eigenvalues corresponds to n = 0.
end