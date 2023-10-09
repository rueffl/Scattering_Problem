function T_matrix = getT_finite(kn, N, lij)
%GETT_FINITE 
%   kn:     wavenumber outside the resonators
%   N:      number of resonators
%   lij:    spacing between resonators

    T_matrix = zeros(2*N);
    T_matrix(1,1) = sqrt(-1)*kn; T_matrix(end,end) = sqrt(-1)*kn;
    if N > 1
        for j = 1:N-1
            T_matrix(2*j:2*j+1,2*j:2*j+1) = makeAl(kn,lij(j));
        end
    end
end

function Al = makeAl(kn,l)
    Al = kn*[-cos(kn*l)/sin(kn*l), 1/sin(kn*l);
          1/sin(kn*l),-cos(kn*l)/sin(kn*l)];
end