function [all_Ln] = get_Lambdas_N(x, k_tr, w, Omega, rs, ks, vr, v0, z, sol, vin_l, vin_r)
%GET_LAMBDAS   Computes the coefficients \Lambda_n, for all n=-K,...,K for
%              a system of N resonators.
%   x:      spatial coordinate
%   k_tr:   truncation parameter
%   w:      quasiwavefrequency
%   Omega:  periodicity of time modulations
%   rs:     Fourier coefficients of \rho_i(t)
%   ks:     Fourier coefficients of \kappa_i(t)
%   vr:     wave speed inside the resonators
%   v0:     wave speed outside the resonators
%   z:      resonator centre
%   sol:    vector of coefficients a_j and b_j
%   vin_l:  left incident wave field, fct of x and n
%   vin_r:  right incident wave field, fct of x and n

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
        if x <= z % left of the resonator
            Ln = Ln - vin_l(z-0.00001,n);
        elseif x >= z % right of the resonator
            Ln = Ln - vin_r(z+0.00001,n);
        end
        all_Ln(n+k_tr+1) = 2*sqrt(-1).*kn.*Ln;
%         all_Ln(n+k_tr+1) = Ln;
    end

end