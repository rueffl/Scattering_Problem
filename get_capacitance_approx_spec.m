%% Numerically compute the Floquet exponents of the NxN system of Capacitance ODEs
%  GCM\Psi + D d/dt \Psi + d/dt 1/\kappa d/dt \Psi = 0. 1/kappa has a
%  finite Fourier series of length 1
%
% Rewrite to the system of 1st order ODEs
% dt/dt \psi_1 = \kappa*psi_2
% dt/dt \psi_2 = -GCM*psi_1 - D/\kappa*psi_2
% Then solve spectrally: d/dt = 1i(\omega + n\Omega) which gives the
% eigenvalue problem
% \omega \psi_1 = - n\Omega \psi_1 - 1i*\kappa*\psi_2
% \omega \psi_2 = - n\Omega \psi_2 - 1i*GCM*\psi_1 - 1i*D*\kappa*psi_2
% 

function w_out = get_capacitance_approx_spec(epsilon_kappa,phase_kappa,Omega,delta,li,v0,vr,C,k_tr)
    GCM = delta*vr^2*diag(1./li)*C;

    M = 1; % Number of Fourier coefficients of 1/\kappa
    N = size(GCM,1);
    if N==1
        d = 2;
    else 
        d = 1;
    end
    
    N_fourier = k_tr; % Length of Fourier series approx
    
    K_mod = zeros(2*M+1,N);
    for i = 1:N
        K_mod(:,i) = [epsilon_kappa/2*exp(-1i*phase_kappa(i)); 1; epsilon_kappa/2*exp(1i*phase_kappa(i))]; % Fourier coefficients of 1/\kappa
    end

    ns = -N_fourier:N_fourier;

    NN = 2*N_fourier+1;
    O = diag(ns.'*Omega);
    e = ones(NN,1);
    INN = eye(NN);
    IN = eye(N);
    iK = zeros(NN*N);
    KD = zeros(NN*N);
    for i = 1:N
        k = zeros(NN,NN);
        for m = -M:M
            k = k+diag(e(1:NN-abs(m))*K_mod(m+M+1,i),m);
        end
        Ii = (i-1)*NN+1:i*NN;
        ik = inv(k);
        iK(Ii,Ii) = ik; %% Fourier coefficients of \kappa
        if (i == 1) || (i == N)
            KD(Ii,Ii) = d*delta*vr^2/(v0*li(i))*ik;
        end
    end

    Z = zeros(NN*N);
    mat = -[kron(IN,O), Z; Z, kron(IN,O)]  - 1i*[Z, iK; -kron(GCM,INN), KD]; % Kroenecker product to get the RHS matrix

    w_out = eig(mat);%,2*N,'smallestabs'); % The eigenvalues of "mat" are approximately \omega + n\Omega for |n| < N_fouier. Taking the smallest eigenvalues corresponds to n = 0.
    w_out = sort(mink(w_out,2*N));
end