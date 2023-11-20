function [N_kn] = getN(kn,N,lij)
%GETN  builds the block-diagonal matrix \mathcal{N}
%   kn:     n-th wave number
%   N:      number of resonators
%   lij:    spacing between the resonators

    N_kn = zeros(2*N);
    N_kn(1,1) = -sqrt(-1)*kn;
    N_kn(end,end) = -sqrt(-1)*kn;
    for j = 1:N-1
        N_kn(2*j:2*j+1,2*j:2*j+1) = makeNl(kn,lij(j));
    end

end

function Al = makeNl(kn,l)
    Al = kn*[cos(kn*l)/sin(kn*l), -1/sin(kn*l);
          1/sin(kn*l),-cos(kn*l)/sin(kn*l)];
end