function [u] = operator_S(x, N, xim, xip, lij, k_tr, k, w, Omega, rs, ks, vr, sol, n, kin, v_in)
%OPERATOR_S operator evaluating the scattered wave at x.
%   x:      evaluation point
%   N:      number of resonators
%   xim:    left boundary points of all resonators
%   xip:    right boundary points of all resonators
%   lij:    spacing between the resonators
%   k_tr:   truncation parameter
%   k:      wavenumber outside the resonators
%   w:      quasiwavefrequency \omega
%   Omega:  period of time modulations
%   rs:     Fourier coefficients of \rho_i(t)
%   ks:     Fourier coefficients of \kappa_i(t)
%   vr:     wave speed inside the resonators
%   sol:    coefficients of the interior solution
%   n:      mode
%   kin:    wave number of incident wave
%   v_in:   n-th mode of the incident wave field


    if x > xip(end)
        i = N;
        a = 0;
        for j = -k_tr:k_tr
            l = ll(j, i, k_tr, w, Omega, rs, ks, vr);
            v = vs(j, n, i, k_tr, w, Omega, rs, ks, vr);
            a = a + ((sol(2*N*(k_tr-j)+2*i-1)*exp(sqrt(-1)*l*xip(end))+sol(2*N*(k_tr-j)+2*i)*exp(-sqrt(-1)*l*xip(end)))*v); 
        end
        a = a - v_in(xip(end)+0.000001,n); % add the influence of the incident wave in the continuity condition
        u = a*exp(sqrt(-1)*k*(x-xip(end)));
    elseif x < xim(1)
        i = 1;
        b = 0;
        for j = -k_tr:k_tr
            l = ll(j, i, k_tr, w, Omega, rs, ks, vr);
            v = vs(j, n, i, k_tr, w, Omega, rs, ks, vr);
            b = b + ((sol(2*N*(k_tr-j)+2*i-1)*exp(sqrt(-1)*l*xim(1))+sol(2*N*(k_tr-j)+2*i)*exp(-sqrt(-1)*l*xim(1)))*v);
        end
        b = b - v_in(xim(1)-0.0000001,n); % add the influence of the incident wave in the continuity condition
        u = b*exp(sqrt(-1)*k*(xim(1)-x));
    else
        for i = 1:(N-1)
            xm = xim(i+1); xp = xip(i);
            if x >= xp && x <= xm
                ab = operator_M(sol, i, k_tr, n, xim, xip, lij, k, w, Omega, rs, ks, vr, v_in);
                u = ab(1)*exp(sqrt(-1)*k*x)+ab(2)*exp(-sqrt(-1)*k*x);
            end
        end
        for i = 1:N
            xm = xim(i); xp = xip(i);
            if x >= xm && x <= xp
                u = 0;
                for j = -k_tr:k_tr
                    l = ll(j, i, k_tr, w, Omega, rs, ks, vr);
                    v = vs(j, n, i, k_tr, w, Omega, rs, ks, vr);
                    u = u + (sol(2*N*(k_tr-j)+2*i-1)*exp(sqrt(-1)*l*x)+sol(2*N*(k_tr-j)+2*i)*exp(-sqrt(-1)*l*x))*v;
                end
            end
        end
    end

end

function lambda = ll(k, i, k_tr, w, Omega, rs, ks, vr)
    C = getC(k_tr, i, w, Omega, rs, ks, vr);
    [~,lambdas] = eig(C,'vector');
    lambdas = sqrt(lambdas);
    lambda = lambdas(k_tr-k+1);
end

function v = vs(k, n, i, k_tr, w, Omega, rs, ks, vr)
    C = getC(k_tr, i, w, Omega, rs, ks, vr);
    [fi,~] = eig(C,'vector');
    if k_tr-n+1 < 1 || k_tr-n+1 > 2*k_tr+1
        v = 0;
    else
        v = fi(k_tr-n+1,k_tr-k+1);
    end
end