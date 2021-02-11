function [C]=compute_plik(Y,ind0,Fct,rc,SUF,Ypos)

% Compute the log likelihood 
% 
% INPUT:
% Y        : Histograms of photon count
% ind0     : set of non-empty time bins
% Fct      : Impulse response functions
% rc       : current reflectivity profile
% SUF      : stored sum for Poisson likelihood computation

% OUTPUT:
% C        : Log likelihood
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen, "Stochastic EM
%       algorithm for fast analysis of single waveform multi-spectral Lidar 
%       data," 2020 28th European Signal Processing Conference (EUSIPCO), 
%       Amsterdam, 2021, pp. 2413-2417, doi: 10.23919/Eusipco47968.2020.9287414



% M vector of size Nx
% SUF matrix of size Nx * Ns
% rc matrix of size Ns * N
% ff matrix of size ii * Nx (Attention ii n'est pas de même taille pour 
%                            tout les n, d'où l'utilisation de cellules matlab)
%                            de taille T * Nx si ii pas dispo
% Fct tensor of size T * Nx * Ns
% Y of size N * ii (ou T si ii n'est pas dispo)
% ind0 cellule de taille N (chaque element d'une cellulle peut etre
%                           n'importe quoi: un vecteur, une matrice, ou
%                           même d'autre cellules! Ici ce sont des vecteurs
%                           de tailles ii: les bin dans lesquels on a des
%                           observations
% Ypos cellule de taille N : chacune de taille ii

% tic
N = size(Y,1);
[~,Nx,Ns]=size(Fct);
C=zeros(N,Nx);
for n=1:N
    ii=ind0{n};
    if ~isempty(ii)
        M = SUF*rc(:,n);
        ff = zeros(length(ii),Nx);
            for ns = 1:Ns
                ff = ff + Fct(ii,:,ns)*rc(ns,n); % to optimize
            end
        C(n,:) = sum(Ypos{n}'.*log(ff)- ones(length(ii),1).*M');
    end
end


% Using mex files inside the loop on N

% C=zeros(N,Nx);
% for n=1:N
%     ii=ind0{n};
%     if ~isempty(ii)
%         M = SUF*rc(:,n);
%         C(n,:) = testii(Ypos{n},rc(:,n),Fct(ii,:,:),M);
%     end
% end


% Using mex file without loop

% C2 = Compute_cplik(Y,rc,Fct,SUF);
