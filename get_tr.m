function [tns,rns] = get_tr(k_tr, w_op, Omega, rs, ks, vr, xim, xip, sol, lij, vin_l, vin_r)
%GET_TR   Computes the transmission and reflection coefficients
%   k_tr:   truncation parameter
%   w_op:   quasiwavefrequency
%   Omega:  periodicity of time modulations
%   rs:     Fourier coefficients of \rho_i(t)
%   ks:     Fourier coefficients of \kappa_i(t)
%   vr:     wave speed inside the resonators
%   xim:    left boundary points of the resonators
%   xip:    right boundary points of the resonators
%   sol:    vector of coefficients a_j and b_j
%   lij:    spacing between the resonators
%   vin_l:  left incident wave field, function of (x,n)
%   vin_r:  right incident wave field, function of (x,n)

    as = sol(1:2:end-1); bs = sol(2:2:end); % coefficients of interior solution
    N = length(xim);

    tns = zeros(2*k_tr+1,N); % transmission coefficients
    rns = zeros(2*k_tr+1,N); % reflection coefficients

    for n = -k_tr:k_tr
        kn = (w_op+n*Omega)/vr;
        for i = 1:N
            if i == 1
                C = getC(k_tr, i, w_op, Omega, rs, ks, vr);
                [f,lambdas] = eig(C,'vector');
                lambdas = sqrt(lambdas);
                fn = f(k_tr-n+1,:);
                beta1 = 0;
                for j = -k_tr:k_tr
                    beta1 = beta1 + (sol(2*N*(k_tr-j)+2*i-1)*exp(sqrt(-1)*lambdas(k_tr+1-j)*xim(1))+sol(2*N*(k_tr-j)+2*i)*exp(-sqrt(-1)*lambdas(k_tr+1-j)*xim(1)))*fn(k_tr+1-j);
                end
                tns(n+k_tr+1,i) = exp(sqrt(-1)*kn*xim(1))*(beta1-vin_l(xim(1)-0.000001,n));
            end
            if i == N
                C = getC(k_tr, i, w_op, Omega, rs, ks, vr);
                [f,lambdas] = eig(C,'vector');
                lambdas = sqrt(lambdas);
                fn = f(k_tr-n+1,:);
                alphaN = 0;
                for j = -k_tr:k_tr
                    alphaN = alphaN + (sol(2*N*(k_tr-j)+2*i-1)*exp(sqrt(-1)*lambdas(k_tr+1-j)*xip(end))+sol(2*N*(k_tr-j)+2*i)*exp(-sqrt(-1)*lambdas(k_tr+1-j)*xip(end)))*fn(k_tr+1-j);
                end
                rns(n+k_tr+1,i) = exp(-sqrt(-1)*kn*xip(end))*(alphaN+vin_r(xip(end)+0.000001,n));
            end
            if i>1 && i<N
                alpha_beta = operator_M(sol, i, k_tr, n, xim, xip, lij, kn, w_op, Omega, rs, ks, vr, @(x,n) 0);
                rns(n+k_tr+1,i) = alpha_beta(1);
                tns(n+k_tr+1,i) = alpha_beta(2);
            end
        end 
    end

end