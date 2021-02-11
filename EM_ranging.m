function [r_out,T0]=EM_ranging(Y,Fc,Ns,Nrow,Ncol,c,si,step_Nx,stepg,nb_label)

% Main algorithm: estimate depth and reflectivity profiles
% 
% INPUT:
% Y         : Histograms of photon count
% Fc        : Impulse response functions
% Ns        : Number of spectral component +1 for the background
% Nrow      : Number of rows
% Ncol      : Number of columns
% c         : Depth TV prior weight
% si        : Neighborhood for patch k-means
% step_Nx   : Likelihood subsampling factor
% stepg     : Gradient step 
% nb_label  : Number of cluster for patch k-means   

% OUTPUT:
% r_out    : Reflectivity estimate
% T0       : Depth estimate
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen, "Stochastic EM
%       algorithm for fast analysis of single waveform multi-spectral Lidar 
%       data," 2020 28th European Signal Processing Conference (EUSIPCO), 
%       Amsterdam, 2021, pp. 2413-2417, doi: 10.23919/Eusipco47968.2020.9287414


%% Initialisation
[N,T]=size(Y);% nb of pixels x histogram length
[~,Nx,~] = size(Fc); % nb of admissible length x nb of spectral component
rc=0.5*ones(Ns,N); % Initialisation reflectivity
r0c = rc; % memory 
r_out=zeros(Ns,N);
iteEM = 0; % Iteration kept for final estimation
c = c.*step_Nx; % Scaling TV prior parameter
ind0=cell(N,1);Ypos=cell(N,1);np=zeros(N,1);
m_compt=1; % Iteration
err0=1e9; % error between iterations
k=2*ones(1,Ns); % gamma prior parameters
theta=0.5*ones(1,Ns); % gamma prior parameters
% gamma prior parameters

cond = 0; % indicates if clusters are defined
kloc = 2*ones(nb_label,Ns); % C-gamma prior parameters
thetaloc = 0.5*ones(nb_label,Ns); % C-gamma prior parameters
% Storing for further computations
for n=1:N
    ind0{n}=find(Y(n,:)>0); % find non-empty time bins
    Ypos{n}=Y(n,ind0{n}); % find photon counts in these bins
    np(n)=numel(ind0{n}); % find number of non-empty bins
end
% Compute a reduced version of the IRF
if step_Nx ~= 1
    Fct = zeros(T,Nx/step_Nx,Ns);
    for i = 1:Nx/step_Nx
        Fct(:,i,:) = sum(Fc(:,(i-1)*step_Nx+1:i*step_Nx,:),2);
    end
else
    Fct = Fc;
end
SUF = squeeze(sum(Fct,1)); % precomputation for Poisson likelihood


%% Initialization T0,P
% Sampling from posterior P using uniform prior
C=compute_plik(Y,ind0,Fct,rc,SUF,Ypos);
C=C-max(C,[],2)*ones(1,Nx/step_Nx);
P=exp(C);
P=P./(sum(P,2)*ones(1,Nx/step_Nx));
T0=denSampling(1:Nx/step_Nx,P);
T0 = T0.*step_Nx;


%% Main iterative process
% while ((m_compt<=15 && derr>1e-3 && err0>1e-5) || m_compt<=15)
while ((err0>1e-2 || m_compt<=10) && m_compt<=10)
    disp(['Iteration : ',num2str(m_compt)])
    %% Expectation / Stochastic Expectation step
    %%TV-depth regularization
    if m_compt>=3 % first iterations without TV
        [C]=compute_plik(Y,ind0,Fct,rc,SUF,Ypos);
        [~,T0]=Estep_Post_sampleT(C,reshape(T0./step_Nx,Nrow,Ncol),c);% using SEM, sampling over T
        T0 = T0.*step_Nx;
    end

    %% Maximization
    if m_compt<3
%       Perform MLE
        rc=Mstep_Wgamma(Y,Ypos,Fc,rc,k,theta,ind0,T0,stepg,SUF,step_Nx);
    else
        if cond ==0 % define cluster if not done yet
            [~,indk]=Patch_cube_multi(reshape(rc(1:Ns-1,:)',Nrow,Ncol,Ns-1),si,nb_label,2); 
            cond = 1;
        end
        % M-step
        [rc,kloc,thetaloc]=Mstep_Cgamma_updpp(Y,Ypos,Fc,rc,kloc,thetaloc,stepg,ind0,T0,indk);
    end

    %% Compute errors
    err0=sum(sum(abs(rc-r0c)))./min(sum(sum(abs(rc))),sum(sum(abs(r0c))));
    r0c=rc;
    
    if m_compt >=10
        r_out=r_out+rc;
        iteEM = iteEM + 1;
    end
    m_compt=m_compt+1;
end
r_out = r_out ./ iteEM;


