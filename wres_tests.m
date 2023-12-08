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


%% Apply vector search method

% tolerance for exiting the vector search method
tol = 1e-10;

% function to iterate with
p = @(ek) find_imagzero(ek,li,Omega,phase_kappa,delta,vr,v0,lij);

eks = linspace(0,4,800); p_eks = zeros(1,length(eks)); k = 1;
for ek = eks
    p_eks(k) = p(ek);
    k = k+1;
end
% figure();
% plot(eks,p_eks,'-','LineWidth',2);

% find minimum
[min_val,idx] = min(p_eks);
min_eks = eks(idx);

%% Determine regions

w_res = @(epsilon_kappa,idx) imag(omega_epsk(epsilon_kappa,li,Omega,phase_kappa,delta,vr,v0,lij,idx));
all_epsk = linspace(0,4,800);
w_epsk = zeros(2*N,length(all_epsk));

% prepare plot
fig = figure();
hold on

change_idx = [];
change_val = [];
val_old = all_epsk(1);
if N == 1

    idx = 1;
    for i = 1:length(all_epsk)
        w_epsk(idx,i) = w_res(all_epsk(i),idx);
        % find change of sign
        if i > 1
            if w_epsk(idx,i)*w_epsk(idx,i-1) < 0
                change_idx = [change_idx,i];
            end
        end
        change_val = all_epsk(change_idx);
    end

    for i = 1:length(change_val)
        val = change_val(i);
        inBetweenx = [val_old.*ones(1,100),fliplr(val.*ones(1,100))];
        inBetweeny = [linspace(min(w_epsk(idx,:))-0.005,max(w_epsk(idx,:))+0.005,100),fliplr(linspace(min(w_epsk(idx,:))-0.005,max(w_epsk(idx,:))+0.005,100))];
        if w_epsk(idx,change_idx(i)-1) < w_epsk(idx,change_idx(i)) 
            str = '#f3846c'; % red
            color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
            fill(inBetweenx,inBetweeny,color,'FaceAlpha',0.3,'LineStyle','none')
        else
            str = '#74ea91'; % green
            color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
            fill(inBetweenx,inBetweeny,color,'FaceAlpha',0.3,'LineStyle','none')
        end
        if abs(w_epsk(idx,change_idx(i)-1) - w_epsk(idx,change_idx(i))) < 10^(-4)
            str = '#f98204'; % orange
            color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
            plot(val.*ones(1,100),linspace(min(w_epsk(idx,:))-0.005,max(w_epsk(idx,:))+0.005,100),'--','Color',color,markersize=8,linewidth=2)
        end
        val_old = val;
    end
    inBetweenx = [val_old.*ones(1,100),fliplr(all_epsk(end).*ones(1,100))];
    inBetweeny = [linspace(min(w_epsk(idx,:))-0.005,max(w_epsk(idx,:))+0.005,100),fliplr(linspace(min(w_epsk(idx,:))-0.005,max(w_epsk(idx,:))+0.005,100))];
    if w_epsk(idx,change_idx(i)-1) < w_epsk(idx,change_idx(i)) 
        str = '#f3846c'; % red
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        fill(inBetweenx,inBetweeny,color,'FaceAlpha',0.3,'LineStyle','none')
    else
        str = '#74ea91'; % green
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        fill(inBetweenx,inBetweeny,color,'FaceAlpha',0.3,'LineStyle','none')
    end

    if abs(w_epsk(idx,change_idx(i)-1) - w_epsk(idx,change_idx(i))) < 10^(-4)
        str = '#f98204'; % orange
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        plot(val.*ones(1,100),linspace(min(w_epsk(idx,:))-0.005,max(w_epsk(idx,:))+0.005,100),'--','Color',color,markersize=8,linewidth=2)
    end
    plot(all_epsk,w_epsk(idx,:),'.k',markersize=8,linewidth=2)
    xlabel('$\varepsilon_{\kappa}$','Interpreter','latex','fontsize',18)
    ylabel('$\omega$','Interpreter','latex','fontsize',18)

else

    for idx = 1:2*N

        subplot(2,N,idx)
        hold on

        change_idx = []; change_val = [];

        for i = 1:length(all_epsk)
            w_epsk(idx,i) = w_res(all_epsk(i),idx);
            % find change of sign
            if i > 1
                if w_epsk(idx,i)*w_epsk(idx,i-1) < 0
                    change_idx = [change_idx,i];
                    change_val = [change_val,all_epsk(i)];
                end
            end
        end

        for i = 1:length(change_val)
            val = change_val(i);
            inBetweenx = [val_old.*ones(1,100),fliplr(val.*ones(1,100))];
            inBetweeny = [linspace(min(w_epsk(idx,:))-0.005,max(w_epsk(idx,:))+0.005,100),fliplr(linspace(min(w_epsk(idx,:))-0.005,max(w_epsk(idx,:))+0.005,100))];
            if w_epsk(idx,change_idx(i)-1) < w_epsk(idx,change_idx(i)) 
                str = '#f3846c'; % red
                color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
                fill(inBetweenx,inBetweeny,color,'FaceAlpha',0.3,'LineStyle','none')
            else
                str = '#74ea91'; % green
                color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
                fill(inBetweenx,inBetweeny,color,'FaceAlpha',0.3,'LineStyle','none')
            end
            if abs(w_epsk(idx,change_idx(i)-1) - w_epsk(idx,change_idx(i))) < 10^(-4)
                str = '#f98204'; % orange
                color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
                plot(val.*ones(1,100),linspace(min(w_epsk(idx,:))-0.005,max(w_epsk(idx,:))+0.005,100),'--','Color',color,markersize=8,linewidth=2)
            end
            val_old = val;
        end
        inBetweenx = [val_old.*ones(1,100),fliplr(all_epsk(end).*ones(1,100))];
        inBetweeny = [linspace(min(w_epsk(idx,:))-0.005,max(w_epsk(idx,:))+0.005,100),fliplr(linspace(min(w_epsk(idx,:))-0.005,max(w_epsk(idx,:))+0.005,100))];
        if w_epsk(change_idx(i)+5) < 0
            str = '#f3846c'; % red
            color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
            fill(inBetweenx,inBetweeny,color,'FaceAlpha',0.3,'LineStyle','none')
        else
            str = '#74ea91'; % green
            color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
            fill(inBetweenx,inBetweeny,color,'FaceAlpha',0.3,'LineStyle','none')
        end
    
        if abs(w_epsk(idx,change_idx(i)-1) - w_epsk(idx,change_idx(i))) < 10^(-4)
            str = '#f98204'; % orange
            color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
            plot(val.*ones(1,100),linspace(min(w_epsk(idx,:))-0.005,max(w_epsk(idx,:))+0.005,100),'--','Color',color,markersize=8,linewidth=2)
        end
        plot(all_epsk,w_epsk(idx,:),'.k',markersize=8,linewidth=2)
        xlabel('$\varepsilon_{\kappa}$','Interpreter','latex','fontsize',18)
        ylabel('$\omega$','Interpreter','latex','fontsize',18)
        
%         for k = 1:2*N
%             plot(all_epsk,w_epsk(k,:),'.k',markersize=8,linewidth=2)
%         end

    end

end





%% Functions

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
        w_res = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0); % subwavelength resonant frequencies
%         w_res = w_res(real(w_res)>=0); % positive subwavelength resonant frequencies
    else
        w_res = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,li,delta,vr,v0); % non-zero subwavelength resonant frequency
    end
    imag_w_res = max(abs(imag(w_res))); % imaginary part of the swl resonant frequencies

end

function [w_res] = omega_epsk(epsilon_kappa,li,Omega,phase_kappa,delta,vr,v0,lij,idx)
% FIND_IMAGZERO function whose root must be found in order to determine setting st imaginary part of resonant frequency is zero
%   epsilon_kappa:  time-modulation amplitude of kappa
%   li:             size of resonators
%   Omega:          modulation frequency of kappa
%   phase_kappa:    phase of time-modulated kappa
%   delta:          contrast parameter
%   vr:             wave speed inside resonators
%   v0:             wave speed in the background medium
%   lij:            spacing between the resonators
%   idx:            index of omega of interest

    N = length(phase_kappa);
    if N > 1
        C = make_capacitance_finite(N,lij); % capacitance matrix
        w_res = get_capacitance_approx_hot(epsilon_kappa,li,Omega,phase_kappa,delta,C,vr,v0); % subwavelength resonant frequencies
%         w_res = w_res(real(w_res)>=0); % positive subwavelength resonant frequencies
    else
        w_res = [get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,li,delta,vr,v0),0]; % non-zero subwavelength resonant frequency
    end
    w_res = w_res(idx);

end





