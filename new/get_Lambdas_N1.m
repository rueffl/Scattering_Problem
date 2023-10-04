function [all_Ln_xl] = get_Lambdas_N1(z,xm,xp,k_tr,w_op,Omega,rs,ks,vr,v0,delta,k_op,k_in,x,i,sol)
%GET_LAMBDAS_N1 Computes the frequency-scattering coefficients of a dilute
%               single-resonator system
%   z:          all midpoints of resonators
%   xm:         left boundary point of the resonator
%   xp:         right boundary point of the resonator
%   k_tr:       truncation parameter
%   w_op:       operating frequency
%   Omega:      frequency of time-modulation
%   rs:         fourier coefficients of 1/rho(t)
%   ks:         fourier coefficients of 1/kappa(t)
%   vr:         wave speed inside the resonators
%   v0:         wave speed outside the resonators
%   delta:      contrast parameter
%   k_op:       operating wave number
%   k_in:       incident wave number
%   x:          spatial location
%   i:          index of the resonator
%   sol:        coefficients of the interior solution

    if isempty(sol)
        MatcalA = getMatcalA(1,[],xm,xp,k_tr,w_op,Omega,rs,ks,vr,delta,v0); % matrix \mathcal{A}
        MatcalF = getF(k_tr, 1, delta, k_op, k_in, xm); % vector \mathcal{F}
        sol = MatcalA\MatcalF; % solve for the interior coefficients, vector \mathbf{w}
    end
    as = sol(1:2:end-1); bs = sol(2:2:end); % coefficients of interior solution
    zi = z(i);
    G = @(k,x) exp(sqrt(-1)*k*abs(x))./(2*sqrt(-1)*k); % Green's function
    
    all_Ln_xl = zeros(1,2*k_tr+1);
    for n = -k_tr:k_tr
        kn = (w_op+n*Omega)/v0;
        Ln = 0;
        Cis = getC(k_tr,1,w_op,Omega,rs,ks,vr); % matrix C
        [fi,l] = eig(Cis,'vector'); % eigenvector and eigenvalue of C
        l = sqrt(l);
        for j = -k_tr:k_tr
            Ln = Ln + (as(k_tr+1+j)*exp(sqrt(-1)*l(k_tr+1-j)*zi)+bs(k_tr+1+j)*exp(-sqrt(-1)*l(k_tr+1-j)*zi))*fi(n+k_tr+1,k_tr+1+j); % define Lambda_n
        end
        if i == 1
            if n == 0
                if x <= xm
                    Ln = Ln - exp(sqrt(-1)*k_in*xm);
                end
            end
        else
            if x <= xm
                Ln_old = get_Lambdas_N1(z,xm,xp,k_tr,w_op,Omega,rs,ks,vr,v0,delta,k_op,k_in,x,i-1,sol);
                Ln = Ln - Ln_old(n+k_tr+1)*G(kn,x-z(i-1));
            end
        end
        all_Ln_xl(n+k_tr+1) = 2*sqrt(-1).*kn.*Ln;
    end

end