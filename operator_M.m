function [phi] = operator_M(sol, i, k_tr, n, xim, xip, lij, k, w, Omega, rs, ks, vr, v_in)
%OPERATOR_M constructs the coefficients of the exterior solution using the 
%     coefficients of the interior solution.
%   sol:    coefficients of the interior solution
%   i:      i-th resonator
%   k_tr:   truncation parameter
%   n:      considered mode
%   xim:    left boundary points of all resonators
%   xip:    right boundary points of all resonators
%   lij:    spacing between the resonators
%   k:      wavenumber outside the resonators
%   w:      quasiwavefrequency \omega
%   Omega:  period of time modulations
%   rs:     Fourier coefficients of \rho_i(t)
%   ks:     Fourier coefficients of \kappa_i(t)
%   vr:     wave speed inside the resonators
%   v_in:   n-th mode of the incident wave

    m = -1/(2*sqrt(-1)*sin(k*lij(i)));
    R = m.*[exp(-sqrt(-1)*k*xim(i+1)), -exp(-sqrt(-1)*k*xip(i)); -exp(sqrt(-1)*k*xim(i+1)), exp(sqrt(-1)*k*xip(i))];
    v = zeros(2,1);
    N = length(sol)/(4*k_tr+2);
    for j = -k_tr:k_tr
        l = ll(j, i, k_tr, w, Omega, rs, ks, vr);
        lp = ll(j, i+1, k_tr, w, Omega, rs, ks, vr);
        vv = vs(j, n, i, k_tr, w, Omega, rs, ks, vr);
        vvp = vs(j, n, i+1, k_tr, w, Omega, rs, ks, vr);
        v(1) = v(1) + (sol(2*N*(k_tr-j)+2*i-1)*exp(sqrt(-1)*l*xip(i))+sol(2*N*(k_tr-j)+2*i)*exp(-sqrt(-1)*l*xip(i)))*vv;
        v(2) = v(2) + (sol(2*N*(k_tr-j)+2*(i+1)-1)*exp(sqrt(-1)*lp*xim(i+1))+sol(2*N*(k_tr-j)+2*(i+1))*exp(-sqrt(-1)*lp*xim(i+1)))*vvp;
    end
    phi = R*v;

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