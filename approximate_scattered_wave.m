function [u] = approximate_scattered_wave(sol,l,k,kr,x,N,zj)
%APPROXIMATE_SCATTERED_WAVE Uses an approximation formula to define the
% scattered wave valid in the static case
%   sol:    coefficients of the interior solution
%   l:      distance between two resonators
%   k:      wave number outside the resonators
%   kr:     wave number inside the resonators
%   x:      evaluation point
%   t:      evaluation time
%   N:      number of resonators
%   zj:     centers of the resonators


    if x >= zj(end)
        aN = sol(2*N-1); bN = sol(2*N);
        a = exp(-sqrt(-1)*k*zj(end))*(aN*exp(sqrt(-1)*kr(end)*zj(end))+bN*exp(-sqrt(-1)*kr(end)*zj(end)));
        u = a*exp(sqrt(-1)*k*x);
    elseif x <= zj(1)
        a1 = sol(1); b1 = sol(2);
        b = exp(sqrt(-1)*k*zj(1))*(a1*exp(sqrt(-1)*kr(1)*zj(1))+b1*exp(-sqrt(-1)*kr(1)*zj(1)));
        u = b*exp(-sqrt(-1)*k*x);
    else
        for i = 1:(N-1)
            xm = zj(i+1); xp = zj(i);
            if x > xp && x < xm
                [alpha,beta] = approximate_alpha_beta(sol,l,i,k,kr(i));
                u = alpha*exp(sqrt(-1)*k*x)+beta*exp(-sqrt(-1)*k*x);
            end
        end
        for i = 1:N
            z = zj(i); 
            if x ==z
                aj = sol(2*i-1); bj = sol(2*i);
                u = aj*exp(sqrt(-1)*kr(i)*x)+bj*exp(-sqrt(-1)*kr(i)*x);
            end
        end
    end

end