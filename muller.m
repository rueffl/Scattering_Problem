function xplus = muller(initial_guess,N,lij,L,xm,xp,k_tr,Omega,rs,ks,vr,delta,v0)
%MULLER Muller's method to find the smallest eigenvalue of the matrix 
%       \matcal{A} as a function of \omega.
    % define three initial guesses
    xnm2 = initial_guess*(1+0.01);
    xnm1 = initial_guess*(1+0.01*exp(1i*2*pi/3));
    % xn = initial_guess*(1+0.01*exp(1i*4*pi/3));
    xn = initial_guess;
    % tolerance for exiting the mullers method
    epsilon = 1e-9;

    % iteration number
    i = 0;
    xplus = xn;

    % function to iterate with
%     p = @(w) minev(getMatcalA(N,lij,xm,xp,k_tr,w,Omega,rs,ks,vr,delta,v0));
    p = @(w) svds(getMatcalA(N,lij,xm,xp,k_tr,w,Omega,rs,ks,vr,delta,v0),1,"smallest");
    
    % implementation of mullers method
    while (abs(p(xn)) > epsilon)

        q = (xn - xnm1)/(xnm1 - xnm2);
        a = q*p(xn) - q*(1+q)*p(xnm1) + q^2*p(xnm2);
        b = (2*q + 1)*p(xn) - (1+q)^2*p(xnm1) + q^2*p(xnm2);
        c = (1 + q)*p(xn);
                
        r = xn - (xn - xnm1)*((2*c)/(b + sqrt(b^2 - 4*a*c)));
        s = xn - (xn - xnm1)*((2*c)/(b - sqrt(b^2 - 4*a*c)));
        
        if(abs(p(r)) < abs(p(s)))
            xplus = r;
            else
            xplus = s;
        end

        xnm2 = xnm1;
        xnm1 = xn;
        xn = xplus;
        
        i = i + 1;
        disp(strcat('Iteration:  ',num2str(i)))
    end
    
end