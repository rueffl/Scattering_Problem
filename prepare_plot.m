function [us, xs] = prepare_plot(k_tr,li,lij,L,N,Omega,delta,phase_kappa,phase_rho,epsilon_kappa,epsilon_rho,xm,xp,vr,v0)
%PREPARE_PLOT Prepares the vector us of the scattered wave evaluated at xs
%   k_tr:           truncation parameter
%   li:             length of each resonator
%   lij:            spacing between the resonators
%   L:              length of the cell
%   N:              number of resonators
%   Omega:          periodicity of time modulations
%   delta:          contrast parameter
%   phase_kappa:    phase of the time-modulation in kappa
%   phase_rho:      phase of the time-modulation in rho
%   epsilon_kappa:  amplitude of the time-modulation in kappa
%   epsilon_rho:    amplitude of the time-modulation in rho
%   xm:             left boundary points of the resonators
%   xp:             right boundary points of the resonators
%   vr:             wave speed inside the resonators
%   v0:             wave speed outside the resonators

    % Fourier coefficients of rhos and kappas
    rs = [];
    ks = [];
    for j = 1:N
        rs_j = [epsilon_rho*exp(-1i*phase_rho(j)),1,epsilon_rho*exp(-1i*phase_rho(j))];
        ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j)),1,epsilon_kappa*exp(-1i*phase_kappa(j))];
        ks = [ks; ks_j];
        rs = [rs; rs_j];
    end

    % Find quasifrequencies
    w_muller = zeros(N,1);
    
    % compute static case
    C = make_capacitance_finite(N,lij);
    w_static = get_capacitance_approx(0,0,li,Omega,phase_rho,phase_kappa,delta,C);
    
    % compute with muller's method 
    for i = 1:N
        initial_guess = w_static(i);
        [w_muller(i),p0] = search_vector_finite(initial_guess,N,lij,xm,xp,k_tr,Omega,rs,ks,vr,delta,v0);
    end
    
    
    % Compute scattered wave
    ls = 100;
    us = zeros(ls,1);
    xs = linspace(xm(1)-1, xp(N)+1,ls);
    % xs = linspace(xm(1)-1, xp(1),ls);
    
    t = 0;
    w = w_muller(1) + 0.0001;
    k = w/v0;
    
    for i = 1:ls
        for j = 1:N

            w0 = w_muller(j);
            k_0 = w0/v0;

            A = getMatcalA(N, lij, xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0);
            F = getF(k_tr, N, delta, k, k_0, xm);
            sol = linsolve(A,F);
            us(i) = us(i) + get_us(xs(i), t, N, xm, xp, lij, k_tr, k, w, Omega, rs, ks, vr, sol, w0);

        end
    end
end