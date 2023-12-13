function [w_out] = get_omega_exact(delta,len,v0,vr,Omega,epsilon_kappa,phi,T)
%GET_OMEGA_EXACT Solves the capacitance approximation for N=1 exactly
%   delta:          contrast parameter
%   len:            length of the resonator
%   v0:             wave speed outside the resonator
%   vr:             wave speed inside the resonator
%   Omega:          modulation frequency of kappa
%   epsilon_kappa:  modulation amplitude of kappa
%   phi:            phase of kappa
%   T:              period of kappa

    % define the integral of the function a(t)
    A = @(t) log(abs(epsilon_kappa*cos(Omega*t+phi)+1))+(4*delta*vr^2)/...
        (Omega*len*sqrt(-epsilon_kappa-1)*sqrt(epsilon_kappa-1)*v0)*...
        atan((sqrt(epsilon_kappa)*tan(0.5*(Omega*t+phi)))/(sqrt(-epsilon_kappa-1)));
    
    % compute omega
    w_out = 1/(sqrt(-1)*T)*mod(A(T)-A(0),2*pi*sqrt(-1));

end