function [F] = getF(k_tr, N, delta, k, k_0, xim)
%GETF Summary of this function goes here
%   k_tr:   truncation parameter
%   N:      number of resonators
%   delta:  contrast parameter
%   k:      wave number outside the resonators 
%   k_0:    wave number of incident wave
%   xim:    left boundary points of all resonators

    Fn = zeros(2*N,1);
    k = 0;
    Fn(1) = sqrt(-1)*delta*(k_0)*exp(sqrt(-1)*k_0*xim(1));

    F = [];
    for i = -k_tr:k_tr
        F = [F;Fn];
    end

end