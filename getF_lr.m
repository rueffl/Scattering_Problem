function [F] = getF_lr(k_tr, N, delta, xim, xip, dx_vin_l, dx_vin_r)
%GETF_LR  The vector \mathcal{F} for the setting with an incident wave from left and right
%   k_tr:       truncation parameter
%   N:          number of resonators
%   delta:      contrast parameter
%   k:          wave number outside the resonators 
%   k_0:        wave number of incident wave
%   xim:        left boundary points of all resonators
%   xip:        right boundary points of all resonators
%   dx_vin_l:   derivative of n-th mode of left incident wave field
%   dx_vin_r:   derivative of n-th mode of right incident wave field

    F = [];
    for n = k_tr:(-1):-k_tr
        Fn = zeros(2*N,1);
        for i = 1:N
            Fn(2*i-1) = -delta*dx_vin_l(xim(i)-0.0000001,n);
            Fn(2*i) = delta*dx_vin_r(xip(i)+0.0000001,n);
        end
        F = [F;Fn];
    end



end