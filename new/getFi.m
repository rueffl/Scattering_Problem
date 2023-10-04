function [F] = getFi(k_tr, delta, wres, wop, Lambdas_xneg, z_old, xim, Omega, v0)
%GETFI Defines the vector of the linear system corresponding to D_i in the
%      dilute regime.
%   k_tr:       truncation parameter
%   N:          number of resonators
%   delta:      contrast parameter
%   wres:       resonant frequency
%   wop:        operating frequency
%   Lambdas:    frequency-scattering coefficients corresponding to D_{i-1}
%   z_old:      centre of the resonator D_{i-1}
%   xim:        left boundary points of all resonators
%   Omega:      frequency of time-modulation
%   v0:         wave speed outside the resonators

    F = zeros(2*(2*k_tr+1),1);
    d = abs(xim-z_old);
    for j = 1:(2*k_tr+1)
        kj = (wop+j*Omega)/v0;
        F(2*j-1) = delta*Lambdas_xneg(j)*exp(sqrt(-1)*kj*d)/2;%(2*(wop-wres));
    end

end