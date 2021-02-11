function rc=Mstep_Wgamma(Y,Ypos,Fc,rc,k,theta,ind0,T0,stepg,SUF,step_Nx)

% Compute the M-step using weak-gamma prior
% 
% INPUT:
% Y           : Histograms of photon count
% Ypos        : photon counts in non-empty time bins
% Fc          : Impulse response functions
% rc          : current reflectivity profile
% k           : gamma prior hyperparameter
% theta       : gamma prior hyperparameter
% ind0        : set of non-empty time bins
% T0          : current depth profile
% stepg       : gradient step
% 
% OUTPUT:
% rc         : reflectivity estimate
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen, "Stochastic EM
%       algorithm for fast analysis of single waveform multi-spectral Lidar 
%       data," 2020 28th European Signal Processing Conference (EUSIPCO), 
%       Amsterdam, 2021, pp. 2413-2417, doi: 10.23919/Eusipco47968.2020.9287414


[N,~]=size(Y);

parfor n=1:N
    ii=ind0{n};
    if isempty(ii)
        rc(:,n)=(k-1)./theta; % mode of gamma distribution - if no photon
                              % detected, max the porterior remains to max
                              % the prior
    else
%     y=Ypos{n};
    FT0 = squeeze(Fc(:,T0(n),:));
%     rc(:,n) = my_NewRa_MLE_beta_multi_m(y,rc(:,n),k,theta,squeeze(Fc(:,T0(n),:)),ii);
    rc(:,n) = test_mle(Y(n,:),rc(:,n),k,theta,FT0,stepg);
%     rc(:,n) = SPIRAL(Y(n,:),rc(:,n),k',theta',squeeze(Fc(:,T0(n),:)),ii,SUF(T0(n)/step_Nx,:));
    end
end


