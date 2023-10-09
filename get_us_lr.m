function [us] = get_us_lr(x, t, N, xim, xip, lij, k_tr, v0, w, Omega, rs, ks, vr, sol, w_res, kin, vin_l, vin_r)
%GET_US_LR Get the scattered wave evaluated at x at time t for the case of left & right incident wave fields.
%   x:      evaluation point
%   t:      evaluation time
%   N:      number of resonators
%   xim:    left boundary points of all resonators
%   xip:    right boundary points of all resonators
%   lij:    spacing between the resonators
%   k_tr:   truncation parameter
%   v0:     wave speed outside of the resonators
%   w:      operating frequency \omega
%   Omega:  period of time modulations
%   rs:     Fourier coefficients of \rho_i(t)
%   ks:     Fourier coefficients of \kappa_i(t)
%   vr:     wave speed inside the resonators
%   sol:    coefficients of the interior solution
%   w_res:  subwavelength resonant frequencies \omega_j, j=1,...,N
%   kin:    wave number of incident wave
%   vin_l:  modes of the left incident wave field
%   vin_r:  modes of the right incident wave field

    us = 0;
    for j = 1:N
        for n = -k_tr:k_tr
            v_in = @(x) vin_l(x,n) + vin_r(x,n);
            kn = (w+n*Omega)/v0;
            us = us + operator_S(x, N, xim, xip, lij, k_tr, kn, w, Omega, rs, ks, vr, sol, n, kin, v_in)*exp(sqrt(-1)*(n*Omega+w_res(j))*t);
        end
    end

end