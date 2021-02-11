function [rMAP,indk]=Patch_cube_multi(rc,si,nb_label,met)

% Compute the patches to initialize the C-gamma prior model
% 
% INPUT:
% rc           : current reflectivity profile
% si           : Neighborhood for patch k-means
% nb_label     : Number of cluster for patch k-means
% met          : for image border extension: 1=symmetric | 2=periodic
% 
% OUTPUT:
% rMAP         : Clustered image, useful for plot checking the clustering
% indk         : Indices of pixel in each patch
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen, "Stochastic EM
%       algorithm for fast analysis of single waveform multi-spectral Lidar 
%       data," 2020 28th European Signal Processing Conference (EUSIPCO), 
%       Amsterdam, 2021, pp. 2413-2417, doi: 10.23919/Eusipco47968.2020.9287414


% Init
% si:  neighborhood of each pixel of size (si*2+1)x(si*2+1)
[Nrow,Ncol,Ns] = size(rc);
N = Nrow*Ncol;
Pat = zeros(Nrow,Ncol,Ns*(2*si+1)^2);

if met == 1 % Symetric
    Pat_temp = [rc(:,si:-1:1,:) rc rc(:,end-1:-1:end-si,:)];
    Pat_temp = [Pat_temp(si:-1:1,:,:); Pat_temp; Pat_temp(end-1:-1:end-si,:,:)];
else % Periodic
    Pat_temp = [rc(:,end-si:end,:) rc rc(:,1:si,:)];
    Pat_temp = [Pat_temp(end-si:end,:,:); Pat_temp; Pat_temp(end-si:end,:,:)];
end
% Init patches(cubes)

for i = si+1:Nrow+si
    for j = si+1:Ncol+si
        tem = Pat_temp(i-si:i+si,j-si:j+si,:);
        Pat(i-si,j-si,:) = tem(:);    
    end
end
Pat = reshape(Pat,N,size(Pat,3));

id = kmeans(Pat,nb_label);
Idx1 = id(:);

% Output treatment
indk = cell(1,nb_label);
rMAP = zeros(Nrow*Ncol,1);
for i = 1:nb_label
    [indk{i}, ~] = find(Idx1==i);
    rMAP(indk{i}) = i;
end

% reshape
% rMAP = reshape(rMAP,sqrt(length(rMAP)),sqrt(length(rMAP)));

