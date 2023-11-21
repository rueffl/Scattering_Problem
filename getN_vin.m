function [N_vin] = getN_vin(k_tr, N, delta, xim, xip, w_op, Omega, v0, lij, vin_l, vin_r)
%GETF_NVIN  The matrix-vector product N x v^{in} which needs to be added to the vector F
%   k_tr:       truncation parameter
%   N:          number of resonators
%   delta:      contrast parameter
%   xim:        left boundary points of all resonators
%   xip:        right boundary points of all resonators
%   w_op:       operating frequency
%   Omega:      frequency of time-modulation
%   v0:         wave speed outside resonators
%   lij:        spacing between resonators
%   vin_l:      left incident wave field, function of (x,n)
%   vin_r:      right incident wave field, function of (x,n)

    N_matrix = [];
    vin = [];
    for n = k_tr:(-1):-k_tr

        kn = (w_op+n*Omega)/v0;
        N_kn = getN(kn,N,lij);
        N_matrix = blkdiag(N_matrix,N_kn);

        vin_n = [];
        for i = 1:N
            vin_n = [vin_n;vin_l(xim(i)-0.00001,n);vin_r(xip(i)+0.00001,n)];
        end
        vin = [vin;vin_n];

    end

    N_vin = delta.*(N_matrix*vin);

end