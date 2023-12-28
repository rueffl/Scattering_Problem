function MatcalA = getMatcalA(N,lij,xm,xp,k_tr,w,Omega,rs,ks,vr,delta,v0)
%GETMATCALA construct the matrix \matcal{A}
%   N:      number of resonators
%   lij:    spacing between the resonators
%   xim:    left boundary points of all resonators
%   xip:    right boundary points of all resonators
%   k_tr:   truncation parameter
%   w:      quasiwavefrequency
%   Omega:  periodicity of time modulations
%   rs:     Fourier coefficients of \rho_i(t)
%   ks:     Fourier coefficients of \kappa_i(t)
%   vr:     wave speed inside the resonators
%   delta:  contrast parameter
%   v:      wave speed outside the resonators

    MatcalA = [];
    %  compute all the Ci and their eigenvectors/values
    fis  = zeros(2*k_tr+1,2*k_tr+1,N);
    list_lambdas = zeros(2*k_tr+1,N);
    for i = 1:N
        Cis = getC(k_tr,i,w,Omega,rs,ks,vr);
        [fi,lambdas_tilde] = eig(Cis,'vector');
        lambdas = sqrt(lambdas_tilde);
        fis(:,:,i) = fi;
        list_lambdas(:,i) = lambdas;
    end
    
    for n = k_tr:(-1):-k_tr
        An = [];
%  define Dirchlet Neumann map
        kn = (w+n*Omega)/v0;
        T_matrix = getT_finite(kn,N,lij);

        for k = 1:k_tr*2+1
%           assemble matrix An (eq35)
            [Gnk,Vnk] = getG_and_V(xm,xp,n,k_tr,rs,fis(:,k,:),list_lambdas(k,:));
            An = [An Gnk-delta*T_matrix*Vnk];
        end
%   assemble matrix MatcalA (eq37)
    MatcalA = [MatcalA; An];
    end

end