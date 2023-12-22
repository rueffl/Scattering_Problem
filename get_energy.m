function [tot_en] = get_energy(tns,rns,k_tr)
%GET_ENERGY Computes the energy of a system using the transmission and reflection coefficient of each mode
%   tns:    transmission coefficient of each mode
%   rns:    reflection coefficients of each mode
%   w_op:   operating wave frequency
%   Omega:  time-modulation frequency
%   v0:     wave speed outside of the resonators
%   xm:     left boundary of the first resonator
%   xp:     right boundary of the last resonator
%   k_tr:   truncation parameter

%     r_sum = 0; t_sum = 0;
%     w_i = imag(w_op);
%     for n = -k_tr:k_tr
%         r_sum = r_sum + rns(n+k_tr+1)*exp(-sqrt(-1)*n*Omega/v0*xm)*exp(sqrt(-1)*n*Omega*t);
%         t_sum = t_sum + tns(n+k_tr+1)*exp(-sqrt(-1)*n*Omega/v0*xp)*exp(sqrt(-1)*n*Omega*t);
%     end
%     tot_en = abs(exp(w_i*t)).^2.*(abs(exp(w_i*xm/v0)).^2.*abs(r_sum).^2+abs(exp(w_i*xp/v0)).^2.*abs(t_sum).^2);

    tot_en = 0;
    for n = -k_tr+k_tr
        tot_en = tot_en + abs(rns(n+k_tr+1)).^2 + abs(tns(n+k_tr+1)).^2;
    end

end