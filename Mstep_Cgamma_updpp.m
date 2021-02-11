function [rc,kloc,thetaloc]=Mstep_Cgamma_updpp(Y,Ypos,Fc,rc,kloc,thetaloc,stepg,ind0,T0,indk)

% Compute the M-step using weak-gamma prior
% 
% INPUT:
% Y           : Histograms of photon count
% Ypos        : Photon counts in non-empty time bins
% Fc          : Impulse response functions
% rc          : Current reflectivity profile
% kloc        : C-gamma prior hyperparameter
% thetaloc    : C-gamma prior hyperparameter
% stepg       : Gradient step
% ind0        : Set of non-empty time bins
% T0          : Current depth profile
% indk        : Indices of pixel in each patch
% 
% OUTPUT:
% rc          : Reflectivity estimate
% kloc        : Updated C-gamma prior hyperparameter
% thetaloc    : Updated C-gamma prior hyperparameter
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen, "Stochastic EM
%       algorithm for fast analysis of single waveform multi-spectral Lidar 
%       data," 2020 28th European Signal Processing Conference (EUSIPCO), 
%       Amsterdam, 2021, pp. 2413-2417, doi: 10.23919/Eusipco47968.2020.9287414


Ns = size(rc,1);

for nbl = 1:length(indk)% clusters
    ind = indk{nbl}; % pixels belonging to current clusters
    for in = 1:length(ind) 
            ii=ind0{ind(in)};
            if isempty(ii)
                % mode of gamma distribution - if no photon detected,
                % max the porterior remains to max the prior
                rc(:,ind(in))=(kloc(nbl,:)-1)/thetaloc(nbl,:); 
            else

            rc(:,ind(in)) = test_mle(Y(ind(in),:),rc(:,ind(in)),kloc(nbl,:),thetaloc(nbl,:),squeeze(Fc(:,T0(ind(in)),:)),stepg);
            end
    end
    [kloc(nbl,1:Ns-1),thetaloc(nbl,1:Ns-1)]=update_Hpprior(rc(1:Ns-1,ind),kloc(nbl,1:Ns-1),thetaloc(nbl,1:Ns-1));
end
