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
    k = -2*delta*vr^2/(v0*len);
    A_nenner = @(t) 2*k*atanh((epsilon_kappa-1)*tan(0.5*(Omega*t+phi))/sqrt(epsilon_kappa^2-1));
    A_denom = @(t) Omega*sqrt(epsilon_kappa^2-1);
    A = @(t) -A_nenner(t)/A_denom(t);

    
    % compute omega
    w_out = 1/(sqrt(-1)*T)*mod(A(T)-A(0),2*pi*sqrt(-1));

end