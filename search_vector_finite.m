function [w0,p0] = search_vector_finite(initial_guess,N,lij,xm,xp,k_tr,Omega,rs,ks,vr,delta,v0)
%SEARCH_VECTOR Searches for the minimum of p in a neighborhood of initial_guess
%   initial_guess : Solution to the static problem
%   alpha :         Quasiperiodicity
%   N :             Number of resonators per cell
%   lij :           Spacing between the resonators
%   L :             Length of the unit cell
%   xm :            Left bounary points of the resonators
%   xp :            Right boundary points of the resonators
%   k_tr :          Truncation parameter
%   Omega :         Frequency of the time-modulations
%   rs :            Fourier coefficients of 1/\rho
%   ks :            Fourier coefficients of 1/\kappa
%   vr :            Wave speed inside the resonators
%   delta :         Contrast parameter 
%   v0 :            Wave speed outside the resonators

    p = @(w) svds(getMatcalA(N,lij,xm,xp,k_tr,w,Omega,rs,ks,vr,delta,v0),1,"smallest");
    steps = 100;
    ws = linspace(initial_guess-2*10^(-4),initial_guess+2*10^(-4),steps);

    ps = zeros(1,steps);
    for i = 1:steps
        ps(i) = p(ws(i));
    end

    [p0,idx] = min(ps);
    w0 = ws(idx);

end