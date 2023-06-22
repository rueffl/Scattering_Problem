function cn = get_cn(sol,kn,fn,lambdas,k_tr,z,n)
%GET_CN Defines the constant c_n
%   sol:        vector of coefficients of the interior solution
%   kn:         wave number outside of the resonator
%   lambdas:    eigenvalues of C_1
%   k_tr:       truncation parameter
%   z:          centre of the resonator
%   n:          n-th resonator is considered

    cn = 0;
    for j = -k_tr:k_tr
        cn = cn + (sol(2*n*(k_tr-j)+1)*exp(sqrt(-1)*lambdas(j+k_tr+1)*z)+sol(2*n*(k_tr-j)+2)*exp(-sqrt(-1)*lambdas(j+k_tr+1)*z))*fn(j+k_tr+1);
    end
    cn = 2*sqrt(-1)*kn*cn;

end