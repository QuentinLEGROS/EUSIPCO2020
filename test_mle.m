function rc = test_mle(y,rc,k,theta,Ft,stepg)

% Perform gradient ascent
% 
% INPUT:
% y           : Histograms of photon count of pixel n
% rc          : current reflectivity profile of pixel n
% k           : gamma prior hyperparameter
% theta       : gamma prior hyperparameter
% Ft          : Impulse response function estimated in the current
%               estimated depth
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


for i=1:1
    %% Main loop
    gt = Ft*rc;% + Bck;
    A = (y'.*Ft)./gt;

    %% Gradient
    %     graf = sum(A - (Ft.*rc'));
    graf = sum(A - Ft,1);

    %% Dirichlet prior
    graf = graf + (k-1)./rc' - (1./theta);

    %% NR update
    rc = rc + stepg*graf';
    rc = min(max(rc,1e-4),1-1e-4);
end
