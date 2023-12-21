function [us] = get_us(x, t, N, xim, xip, lij, k_tr, v0, w_op, Omega, rs, ks, vr, sol, w_res, kin, v_in)
%GET_US Get the scattered wave evaluated at x at time t.
%   x:      evaluation point
%   t:      evaluation time
%   N:      number of resonators
%   xim:    left boundary points of all resonators
%   xip:    right boundary points of all resonators
%   lij:    spacing between the resonators
%   k_tr:   truncation parameter
%   v0:     wave speed outside of the resonators
%   w_op:   operating frequency \omega
%   Omega:  period of time modulations
%   rs:     Fourier coefficients of \rho_i(t)
%   ks:     Fourier coefficients of \kappa_i(t)
%   vr:     wave speed inside the resonators
%   sol:    coefficients of the interior solution
%   w_res:  subwavelength resonant frequencies \omega_j, j=1,...,N
%   kin:    wave number of incident wave
%   v_in:   incident wave field, function of x and n

    us = 0;
    for n = -k_tr:k_tr
        kn = (w_op+n*Omega)/v0;
        us = us + operator_S(x, N, xim, xip, lij, k_tr, kn, w_op, Omega, rs, ks, vr, sol, n, kin, v_in)*exp(sqrt(-1)*(n*Omega-w_op)*t);
    end

end