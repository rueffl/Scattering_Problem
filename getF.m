function [F] = getF(k_tr, N, delta, k, k_0, xim)
%GETF Summary of this function goes here
%   k_tr:   truncation parameter
%   N:      number of resonators
%   delta:  contrast parameter
%   k:      wave number outside the resonators 
%   k_0:    wave number of incident wave
%   xim:    left boundary points of all resonators

    Fn = zeros(2*N,1);
    Fn(1) = sqrt(-1)*delta*(k-k_0)*exp(sqrt(-1)*k_0*xim(1));
    
    F = Fn;
    for i = 1:(2*k_tr)
        F = vertcat(F,Fn);
    end

end