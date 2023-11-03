function [alphas,betas] = get_alpha_beta(sol,lambdas,fs,xim,xip,w_op,Omega,v0,k_tr,v_in_r,v_in_l)
%GET_ALPHA_BETA Computes the coefficients alpha_n^i and beta_n^i for a resonator D_i, for all n
%   sol:        vector with coefficients a_n^i and b_n^i
%   lambdas:    sqrt of the eigenvalues of the matrix C_i
%   fs:         eigenvectors of the matrix C_i
%   xim:        left boundary point of D_i
%   xip:        right boundary point of D_i
%   w_op:       operating wave frequency
%   Omega:      frequency of time-modulation
%   v0:         wave speed outside D
%   k_tr:       truncation parameter
%   v_in_r:     right incident wave on D_i
%   v_in_l:     left incident wave on D_i

    alphas = zeros(1,2*k_tr+1); betas = zeros(1,2*k_tr+1);

    for n = -k_tr:k_tr

        kn = (w_op+n*Omega)/v0;

        % calculate alpha_n^i and beta_n^i
        alpha = 0; beta = 0;
        for j = -k_tr:k_tr
            l = lambdas(k_tr+1-j);
            v = fs(k_tr+1-n,k_tr+1-j);
            alpha = alpha + ((sol(2*(k_tr-j)+1)*exp(sqrt(-1)*l*xip)+sol(2*(k_tr-j)+2)*exp(-sqrt(-1)*l*xip))*v);
            beta = beta + ((sol(2*(k_tr-j)+1)*exp(sqrt(-1)*l*xim)+sol(2*(k_tr-j)+2)*exp(-sqrt(-1)*l*xim))*v);
        end
        alpha = alpha - v_in_r(xip+0.000001,n); % add the influence of the incident wave in the continuity condition
        alpha = alpha*exp(-sqrt(-1)*kn*xip);
        beta = beta - v_in_l(xim-0.0000001,n); % add the influence of the incident wave in the continuity condition
        beta = beta*exp(sqrt(-1)*kn*xim);

        alphas(n+k_tr+1) = alpha;
        betas(n+k_tr+1) = beta;

    end
end