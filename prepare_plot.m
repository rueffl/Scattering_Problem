function [us, xs] = prepare_plot(k_tr,li,lij,L,N,Omega,delta,phase_kappa,phase_rho,epsilon_kappa,epsilon_rho,xm,xp,vr,v0,xs)
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
    C = make_capacitance_finite(N,lij);
    w_out = get_capacitance_approx(epsilon_kappa,epsilon_rho,li,Omega,phase_rho,phase_kappa,delta,C);
    w_muller = w_out(real(w_out)>=0);
    
    % compute with muller's method 
%     for i = 1:N
%         initial_guess = w_res(i);
%         w_muller(i) = muller(initial_guess,N,lij,xm,xp,k_tr,Omega,rs,ks,vr,delta,v0);
%     end
    
    
    % Compute scattered wave
    if length(xs) == 0
        ls = 100;
        us = zeros(ls,1);
        xs = linspace(xm(1)-1, xp(N)+1,ls);
    else
        ls = length(xs);
        us = zeros(ls,1);
    end
    
    t = 0;
    w = mean(w_muller) + 0.0001;
    k = w_muller(1)/v0;
    
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






% function [us, xs] = prepare_plot(k_tr,li,lij,L,N,Omega,delta,phase_kappa,phase_rho,epsilon_kappa,epsilon_rho,xm,xp,vr,v0)
% %PREPARE_PLOT Prepares the vector us of the scattered wave evaluated at xs
% %   k_tr:           truncation parameter
% %   li:             length of each resonator
% %   lij:            spacing between the resonators
% %   L:              length of the cell
% %   N:              number of resonators
% %   Omega:          periodicity of time modulations
% %   delta:          contrast parameter
% %   phase_kappa:    phase of the time-modulation in kappa
% %   phase_rho:      phase of the time-modulation in rho
% %   epsilon_kappa:  amplitude of the time-modulation in kappa
% %   epsilon_rho:    amplitude of the time-modulation in rho
% %   xm:             left boundary points of the resonators
% %   xp:             right boundary points of the resonators
% %   vr:             wave speed inside the resonators
% %   v0:             wave speed outside the resonators
% 
%     % Fourier coefficients of rhos and kappas
%     rs = [];
%     ks = [];
%     for j = 1:N
%         rs_j = [epsilon_rho*exp(-1i*phase_rho(j)),1,epsilon_rho*exp(-1i*phase_rho(j))];
%         ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j)),1,epsilon_kappa*exp(-1i*phase_kappa(j))];
%         ks = [ks; ks_j];
%         rs = [rs; rs_j];
%     end
% 
%     % Find quasifrequencies
%     w_muller = zeros(2*N,1);
%     
%     % compute static case
%     C = make_capacitance_finite(N,lij);
%     w_static = get_capacitance_approx(0,0,li,Omega,phase_rho,phase_kappa,delta,C);
%     w_static = get_capacitance_approx_spec(0,zeros(1,N),Omega,delta,li,C);
%     w_static = w_static(real(w_static)>=0);
%     w_out = get_capacitance_approx_spec(epsilon_kappa,phase_kappa,Omega,delta,li,C);
%     w_out = w_out(real(w_out)>=0);
%     
%     % compute with muller's method 
%     for i = 1:N
%         initial_guess = w_static(i);
%         [w_muller(i)] = muller(initial_guess,N,lij,xm,xp,k_tr,Omega,rs,ks,vr,delta,v0);
%     end
%     
%     w_muller = w_muller(real(w_muller)>=0);
% 
% %     w_out = get_capacitance_approx_spec(epsilon_kappa,phase_kappa,Omega,delta,li,C);
% %     w_out = w_out(real(w_out)>=0);
%     
%     
%     % Compute scattered wave
%     ls = 100;
%     us = zeros(ls,1);
%     xs = linspace(xm(1)-1, xp(N)+1,ls);
%     % xs = linspace(xm(1)-1, xp(1),ls);
%     
%     t = 0;
%     w = w_muller(1) + 0.0001;
%     k = w/v0;
%     
%     for i = 1:ls
%         for j = 1:N
% 
%             w0 = w_muller(j);
%             k_0 = w0/v0;
% 
%             A = getMatcalA(N, lij, xm, xp, k_tr, w, Omega, rs, ks, vr, delta, v0);
%             F = getF(k_tr, N, delta, k, k_0, xm);
%             sol = linsolve(A,F);
%             us(i) = us(i) + get_us(xs(i), t, N, xm, xp, lij, k_tr, k, w, Omega, rs, ks, vr, sol, w0);
% 
%         end
%     end
% end