function [all_Ln] = get_Lambdas(xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0, k, k_0, z, sol)
%GET_LAMBDAS   Computes the coefficients \Lambda_n, for all n=-K,...,K for
%              the single-resonator case.
%   xim:    left boundary points of all resonators
%   xip:    right boundary points of all resonators
%   k_tr:   truncation parameter
%   w:      quasiwavefrequency
%   Omega:  periodicity of time modulations
%   rs:     Fourier coefficients of \rho_i(t)
%   ks:     Fourier coefficients of \kappa_i(t)
%   vr:     wave speed inside the resonators
%   delta:  contrast parameter
%   v:      wave speed outside the resonators
%   k:      wave number outside the resonators 
%   k_0:    wave number of incident wave
%   z:      resonator centre

    % Compute solution coefficients
    if isempty(sol)
        N = 1; lij = [];
        A = getMatcalA(N, lij, xm(1), xp(1), k_tr, w, Omega, rs, ks, vr, delta, v0); % matrix \mathcal{A}
        F = getF(k_tr, N, delta, k, k_0, xm); % vector \mathcal{F}
        sol = linsolve(A,F); % solve \mathcal{A}\mathbf{w}=\mathcal{F}
    end
    as = sol(1:2:end-1); bs = sol(2:2:end); % coefficients of interior solution
    
    % Compute Lambda_n
    all_Ln = zeros(1,2*k_tr+1);
    for n = -k_tr:k_tr
        kn = (w+n*Omega)/v0;
        Ln = 0;
        Cis = getC(k_tr,1,w,Omega,rs,ks,vr); % matrix C
        [fi,l] = eig(Cis,'vector'); % eigenvector and eigenvalue of C
        l = sqrt(l);
        for j = -k_tr:k_tr
            Ln = Ln + (as(k_tr+1+j)*exp(sqrt(-1)*l(k_tr+1-j)*z)+bs(k_tr+1+j)*exp(-sqrt(-1)*l(k_tr+1-j)*z))*fi(n+k_tr+1,k_tr+1+j); % define Lambda_n
        end
        all_Ln(n+k_tr+1) = 2*sqrt(-1).*kn.*Ln;
    end

end