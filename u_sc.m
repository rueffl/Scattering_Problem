function [usc] = u_sc(x,t,rs,ks,vr,delta,k_op,kin,w_op,Omega,v0,win,k_tr_n,k_tr_m,z,wres)
%U_SC evaluates the scattered wave at x at time t for the dilute regime
%   x:          spatial coordinate
%   t:          temporal coordinate
%   rs:         fourier coefficients of 1/\rho
%   ks:         fourier coefficients of 1/\kappa
%   vr:         wave speed inside the resonators
%   delta:      contrast parameter
%   k_op:       operating wave number
%   kin:        incident wave number
%   w_op:       operating frequency
%   Omega:      modulation frequency
%   v0:         wave speed outside the resonators
%   win:        incident frequency
%   k_tr:       truncation parameter
%   z:          resonators
%   wres:       resonating frequency

    N = length(z);
    u = 0;
    G = @(k,x) exp(sqrt(-1).*k.*abs(x))./(2.*sqrt(-1).*k); % Green's function 
    uin = @(x,t) exp(kin.*x+win.*t); % incident wave

    for i = 1:N
        zi = z(i); % resonator centre
        all_Ln = get_Lambdas(zi-0.000250000000001, zi+0.00025, k_tr_n, w_op, Omega, rs, ks, vr, delta, v0, k_op, kin, zi, []); % all \Lambda_n(z_i)
        for n = -k_tr_n:k_tr_n
            kn = (w_op+n*Omega)/v0;
%             u = u + all_Ln(n+k_tr_n+1)*G(kn,x-zi)*uin(0,t)*exp(sqrt(-1)*n*Omega*t);
            u = u + all_Ln(n+k_tr_n+1)*G(kn,x-zi)*exp(sqrt(-1)*(wres+n*Omega)*t);
%             u = u + all_Ln(n+k_tr_n+1)*G(kn,x-zi)*uin(zi,t)*exp(sqrt(-1)*n*Omega*t); %*exp(sqrt(-1)*(wres+n*Omega)*t);
        end
    end

    usc = u;%./(w-wres);

end