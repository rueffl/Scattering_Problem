function [all_Ln] = get_Lambdas_mn(xm, xp, k_tr_n, k_tr_m, w, Omega, rs, ks, vr, delta, v0, k_0)
%GET_LAMBDAS   Computes the coefficients \Lambda_n, for all n=-K,...,K for
%              the single-resonator case.
%   xim:    left boundary points of all resonators
%   xip:    right boundary points of all resonators
%   k_tr_n: truncation parameter wrt n
%   k_tr_m: truncation parameter wrt m
%   w:      quasiwavefrequency
%   Omega:  periodicity of time modulations
%   rs:     Fourier coefficients of \rho_i(t)
%   ks:     Fourier coefficients of \kappa_i(t)
%   vr:     wave speed inside the resonators
%   delta:  contrast parameter
%   v0:      wave speed outside the resonators
%   k:      wave number outside the resonators 
%   k_0:    wave number of incident wave
    
    % Compute Lambda_{mn}
    N = 1; lij = [];
    all_Ln = zeros(2*k_tr_m+1,2*k_tr_n+1);
    Cis = getC(k_tr_n,1,w,Omega,rs,ks,vr); % matrix C
    [fi,l] = eig(Cis,'vector'); % eigenvector and eigenvalue of C
    for m = -k_tr_m:k_tr_m
        km = (w+m*Omega)/v0;
        A = getMatcalA(N, lij, xm(1), xp(1), k_tr_n, w+m*Omega, Omega, rs, ks, vr, delta, v0); % matrix \mathcal{A}
        F = getF(k_tr_n, N, delta, km, k_0, xm); % vector \mathcal{F}
        sol = linsolve(A,F); % solve \mathcal{A}\mathbf{w}=\mathcal{F}
        as = sol(1:2:end-1); bs = sol(2:2:end); % coefficients of interior solution
        for n = -k_tr_n:k_tr_n
            kn = (w+n*Omega)/v0;
            Ln = 0;
            for j = 1:(2*k_tr_n+1)
                Ln = Ln + (as(j)+bs(j))*fi(n+k_tr_n+1,j); % define Lambda_n
            end
            all_Ln(m+k_tr_m+1,n+k_tr_n+1) = 2*sqrt(-1).*kn.*Ln;
        end
    end

end