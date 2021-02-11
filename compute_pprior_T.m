function pprior=compute_pprior_T(T,Nz,V,c)


% Compute the log of the depth TV prior
% 
% INPUT:
% T        : Current depth estimate
% Nz       : Number of admissible position
% V        : indices of pixel neighbors
% c        : TV prior weight
%
% OUTPUT:
% pprior   : Log TV prior
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen, "Stochastic EM
%       algorithm for fast analysis of single waveform multi-spectral Lidar 
%       data," 2020 28th European Signal Processing Conference (EUSIPCO), 
%       Amsterdam, 2021, pp. 2413-2417, doi: 10.23919/Eusipco47968.2020.9287414


[N,NN]=size(V);
Xm=200;
T1=ones(N,1)*(1:Nz);
pprior=zeros(N,Nz);
for t=1:NN
X=abs(T1-T(V(:,t))*ones(1,Nz));
X(isnan(X))=0;
X(X>Xm)=Xm;
pprior=pprior + X;
end

pprior=-2*c*pprior;