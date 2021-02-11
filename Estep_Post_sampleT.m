function [P,T_out]=Estep_Post_sampleT(plik0,T,c)

% Compute the TV prior and posterior distribution
% 
% INPUT:
% plik0     : Log likelihood
% T         : Current depth estimate
% c         : TV prior weight
%
% OUTPUT:
% P         : Posterior distribution
% T_out     : New depth estimates
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen, "Stochastic EM
%       algorithm for fast analysis of single waveform multi-spectral Lidar 
%       data," 2020 28th European Signal Processing Conference (EUSIPCO), 
%       Amsterdam, 2021, pp. 2413-2417, doi: 10.23919/Eusipco47968.2020.9287414


Nx=size(plik0,2);
[Nrow,Ncol]=size(T);
Z=NaN*ones(Nrow+2,Ncol+2);
P1=zeros(Nrow+2,Ncol+2,Nx);
P1(2:end-1,2:end-1,:)=reshape(plik0,Nrow,Ncol,Nx);
P1=reshape(P1,(Nrow+2)*(Ncol+2),Nx);
T1=NaN*ones(Nrow+2,Ncol+2);
T1(2:end-1,2:end-1)=T;
% construction damier avec 0 et 1 (NAN sur les bords)
Z(2:2:end-1,2:2:end-1)=0;
Z(3:2:end-1,3:2:end-1)=0;
Z(2:2:end-1,3:2:end-1)=1;
Z(3:2:end-1,2:2:end-1)=1;

C=zeros((Nrow+2)*(Ncol+2),Nx);



T2=T1;
ind0=cell(2,1);
for t=1:2
for r=randperm(2)
% vecteur d'indice du damier tq find(damier == 0 ou 1)
ind0{r}=find(Z==r-1);
% voisinage pour chaque point du vecteur d'indice
V=[ind0{r}-1 ind0{r}+1 ind0{r}+Nrow+2  ind0{r}-(Nrow+2)];
% valeur de plik sur le damier
plik=P1(ind0{r},:);
pprior=compute_pprior_T(T2(:),Nx,V,c);
ppost=plik+pprior;
ppost=ppost-max(ppost,[],2)*ones(1,Nx);
post=exp(ppost);
ppost=post./(sum(post,2)*ones(1,Nx));
T2(ind0{r})=denSampling(1:Nx,ppost);
end
end

T_out=T2;
for r=randperm(2)
ind0{r}=find(Z==r-1);
V=[ind0{r}-1 ind0{r}+1 ind0{r}+Nrow+2  ind0{r}-(Nrow+2)];
plik=P1(ind0{r},:);
pprior=compute_pprior_T(T2(:),Nx,V,c);
C((ind0{r}),:)=plik+pprior;
% [~,T_out(ind0{r})]=max(ppost,[],2);
end

C=reshape(C,Nrow+2,Ncol+2,Nx);
C=C(2:end-1,2:end-1,:);
C=reshape(C,Nrow*Ncol,Nx);
C=C-max(C,[],2)*ones(1,Nx);
P=exp(C);
P=P./(sum(P,2)*ones(1,Nx));
T_out=T_out(2:end-1,2:end-1);
T_out=T_out(:);