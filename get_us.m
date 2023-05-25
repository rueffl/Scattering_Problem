function [us] = get_us(x, t, N, xim, xip, lij, k_tr, k, w, Omega, rs, ks, vr, sol, w0)
%GET_US Get the scattered wave evaluated at x at time t.
%   x:      evaluation point
%   t:      evaluation time
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
%   w0:     quasifrequency of incident wave

    us = 0;
    for n = -k_tr:k_tr
        us = us + operator_S(x, N, xim, xip, lij, k_tr, k, w, Omega, rs, ks, vr, sol, n)*exp(sqrt(-1)*(n*Omega+w0)*t);
    end

end