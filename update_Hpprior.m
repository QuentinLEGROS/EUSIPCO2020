function [k,theta]=update_Hpprior(rc,k,theta)

% Compute the M-step using weak-gamma prior
% 
% INPUT:
% rc        : Current reflectivity profile
% k         : C-gamma prior hyperparameter
% theta     : C-gamma prior hyperparameter
% 
% OUTPUT:
% k         : Updated C-gamma prior hyperparameter for current patch 
% theta     : Updated C-gamma prior hyperparameter for current patch 
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen, "Stochastic EM
%       algorithm for fast analysis of single waveform multi-spectral Lidar 
%       data," 2020 28th European Signal Processing Conference (EUSIPCO), 
%       Amsterdam, 2021, pp. 2413-2417, doi: 10.23919/Eusipco47968.2020.9287414


% gamma hyper prior parameters (for k)
g=1.01;h=0.5;
% inv gamma hyper prior parameters (for theta)
alpha=1.01;beta=0.5;

N=size(rc,2);
slog=sum(log(rc),2);
sn = sum(rc,2);

for rep = 1:50
    for t=1:size(rc,1)
        
        err=10;
        k0=k(t);
        it=1;
        while (err>1e-4 && it <20)
            g = slog(t) - N*psi(k(t)) - N*log(theta(t)) + (g-1)/k(t) - (1/h) ;
            gp = - N*psi(1,k(t))-((g-1)/(k(t)^2));
            k(t)=min(max(k(t)-g/gp,1.01),100);
            err=abs(k(t)-k0);
            k0=k(t);
            it = it+1;
        end
        
        theta0=theta(t);
        err=10;
        it=1;
        while (err>1e-4 && it <20)
            g = sn(t)/(theta(t)^2) - (N*k(t))/theta(t) - (alpha+1)/theta(t) + beta/(theta(t)^2);
            gp = -2*(sn(t)/(theta(t)^3)) + (N*k(t))/(theta(t)^2) + (alpha+1)/(theta(t)^2) - 2*beta/(theta(t)^3);
            theta(t)=min(max(theta(t)-g/gp,0.01),100);
            err=abs(theta(t)-theta0);
            theta0=theta(t);
            it = it+1;
        end
        
    end
end
