function w_out = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,l,delta,vr,v0)
%GET_CAPACITANCE_APPROX_SPEC_IM_N1_1D Computes the subwavelength
%quasifrequencies of a single-resonator system in 1D
%   epsilon_kappa:  modulation amplitude of kappa
%   Omega:          frequency of kappa
%   l:              length of the resonator
%   delta:          contrast parameter
%   vr:             wave speed inside the resonator
%   v0:             wave speed outside the resonator

    M = 1; % Number of Fourier coefficients of 1/\kappa
    N_fourier = 10; % Length of Fourier series approx
    
    K_mod = [epsilon_kappa/2; 1; epsilon_kappa/2]; % Fourier coefficients of 1/\kappa

    ns = -N_fourier:N_fourier;

    NN = 2*N_fourier+1;
    O = diag(ns.'*Omega);
    e = ones(NN,1);
    K = zeros(NN);
    for m = -M:M
        K = K+diag(e(1:NN-abs(m))*K_mod(m+M+1),m);
    end
    iK = inv(K); %% Fourier coefficients of \kappa

    c = 2*delta*(vr)^2/(v0*l);
    mat = -O-1i*c.*iK; 

    %w_out = eigs(mat,1,'smallestreal'); % The eigenvalues of "mat" are approximately \omega + n\Omega for |n| < N_fouier. Taking the smallest eigenvalues corresponds to n = 0.
    w_out = eig(mat);
    [x, I] = min(abs(real(w_out)));
    w_out = w_out(I);
end