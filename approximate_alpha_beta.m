function [alpha,beta] = approximate_alpha_beta(sol,l,j,k,kr)
%APPROXIMATE_ALPHA_BETA approximates alpha and beta in terms of the
%interior coefficients, valid for the static case
%   sol:    coefficients of the interior solution
%   l:      distance between two resonators
%   j:      j-th resonator
%   k:      wave number outside the resonators
%   kr:     wave number inside the resonators

    aj = sol(2*j-1); bj = sol(2*j);
    ajp = sol(2*(j+1)-1); bjp = sol(2*(j+1));

    alpha = 0.5*aj*(-j+1/sqrt(-1)*(j-1)*j*kr*l-1/(2*sqrt(-1))*j^2*k*l+1/(sqrt(-1)*k*l))...
        -0.5*bj*(j+1/sqrt(-1)*(j-1)*j*kr*l+1/(2*sqrt(-1))*j^2*k*l-1/(sqrt(-1)*k*l))...
        +0.5*ajp*(-(j-1)+1/sqrt(-1)*(j-1)*j*kr*l+1/(2*sqrt(-1))*(j-1)^2*k*l-1/(sqrt(-1)*k*l))...
        -0.5*bjp*((j-1)+1/sqrt(-1)*(j-1)*j*kr*l-1/(2*sqrt(-1))*(j-1)^2*k*l+1/(sqrt(-1)*k*l));
    beta = 0.5*aj*(-j+1/sqrt(-1)*(j-1)*j*kr*l+1/(2*sqrt(-1))*j^2*k*l-1/(sqrt(-1)*k*l))...
        -0.5*bj*(j+1/sqrt(-1)*(j-1)*j*kr*l-1/(2*sqrt(-1))*j^2*k*l+1/(sqrt(-1)*k*l))...
        +0.5*ajp*((j-1)-1/sqrt(-1)*(j-1)*j*kr*l-1/(2*sqrt(-1))*(j-1)^2*k*l+1/(sqrt(-1)*k*l))...
        -0.5*bjp*(-(j-1)-1/sqrt(-1)*(j-1)*j*kr*l+1/(2*sqrt(-1))*(j-1)^2*k*l-1/(sqrt(-1)*k*l));

end