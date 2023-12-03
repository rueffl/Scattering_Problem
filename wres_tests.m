%% Set structure setting
clear 
format long

% Settings for the material's structure
k_tr = 4; % truncation parameters as in remark 3.3
N = 2; % number of the resonator
% L = 2000; % length of the domain \mathcal{U}
% spacing = L/N; lij = ones(1,N-1).*spacing; % spacing between the resonators
spacing = 10; lij = ones(1,N-1).*spacing; % spacing between the resonators
len = 2; li = ones(1,N).*len; % length of the resonator
L = sum(li)+sum(lij); % length of the unit cell
Ls = zeros(2*N-1,1);
Ls(1:2:end) = li;
Ls(2:2:end) = lij;
xipm = [0,cumsum(Ls)']-(len*(N/2)+spacing*(N-1)/2); % all boundary points, make sure the resonators are aligned symmetrically wrt 0
xm = xipm(1:2:end); % LHS boundary points
xp = xipm(2:2:end); % RHS boundary points
z = (xm+xp)./2; % centers of resonators
delta = 0.0001; % small contrast parameter
t = 0; % time

vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators

% Setting for the material's time-modulation
Omega = 0.03; % modulation frequency
T = 2*pi/Omega;
phase_kappa = zeros(1,N); % modulation phases of kappa
phase_rho = zeros(1,N); % modulation phases of rho
for i = 1:(N-1)
    phase_kappa(i+1) = pi/i;
    phase_rho(i+1) = pi/i;
end
epsilon_rho = 0; % modulation amplitude of rho


%% Apply Muller's method

% define three initial guesses
initial_guess = 0.1;
xnm2 = initial_guess*(1+0.01);
xnm1 = initial_guess*(1-0.01);
xn = initial_guess;
% tolerance for exiting the mullers method
tol = 1e-10;

% iteration number
i = 0;
xplus = xn;

% function to iterate with
p = @(ek) find_imagzero(ek,li,Omega,phase_kappa,delta,vr,v0,lij);

% eks = linspace(0,4,400); p_eks = zeros(1,length(eks)); k = 1;
% for ek = eks
%     p_eks(k) = p(ek);
%     k = k+1;
% end
% plot(eks,p_eks,'-','LineWidth',2);

% implementation of mullers method
while (abs(p(xn)) > tol)

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
    disp(strcat('Iteration  ',num2str(i),' xn=',num2str(xn),' error=',num2str(p(xn))))
end


%% Function which we want to minimize

function [imag_w_res] = find_imagzero(epsilon_kappa,li,Omega,phase_kappa,delta,vr,v0,lij)
% FIND_IMAGZERO function whose root must be found in order to determine setting st imaginary part of resonant frequency is zero
%   epsilon_kappa:  time-modulation amplitude of kappa
%   li:             size of resonators
%   Omega:          modulation frequency of kappa
%   phase_kappa:    phase of time-modulated kappa
%   delta:          contrast parameter
%   vr:             wave speed inside resonators
%   v0:             wave speed in the background medium
%   lij:            spacing between the resonators

    N = length(phase_kappa);
    if N > 1
        C = make_capacitance_finite(N,lij); % capacitance matrix
        w_muller = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0); % subwavelength resonant frequencies
        w_res = w_muller(real(w_muller)>=0); % positive subwavelength resonant frequencies
    else
        w_res = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,li,delta,vr,v0); % non-zero subwavelength resonant frequency
    end
    imag_w_res = max(abs(imag(w_res))); % imaginary part of the swl resonant frequencies

end






