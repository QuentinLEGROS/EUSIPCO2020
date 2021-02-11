% function [Y,Fc,Nrow,Ncol,A0,T0]=initdata()

B=[0.01 0.05 0.1 0.5 1 5]; % background levels
alp1=[10 25 100 200]; % gain/relates to the overall number of photons

B = B(1);
alp1 = alp1(3);

str = load('test_data.mat');
T0 = str.T0;
A0 = str.A0;

% load init_lego.mat %load ground truth A0 (reflectivity) and T0 (range)
F = load('F_lego_473_640.mat'); % load impulse response 1500*300
F = F.F;

T=size(F,1); % length of the histograms
Fc(:,:,1) = squeeze((F(:,:,1)+eps)./(ones(T,1)*sum(squeeze(F(:,:,1))+eps,1))); % normalised distributions of the signal photons
Fc(:,:,2) = squeeze((F(:,:,2)+eps)./(ones(T,1)*sum(squeeze(F(:,:,2))+eps,1)));
Fc(:,:,3) = squeeze((F(:,:,3)+eps)./(ones(T,1)*sum(squeeze(F(:,:,3))+eps,1)));

Fc(:,:,4) = 1/T; % Background photons p.d.f.
Nbspectr = size(A0,3); % Number of spectral component

Nrow=size(T0,1);
Ncol=size(T0,2);

H1 = A0(:,:,1); 
H2 = A0(:,:,2); 
H3 = A0(:,:,3); 

T0 = ceil(T0);

X0=   alp1*(Fc(:,T0(:),1).*(ones(T,1)*H1(:)')) ...
    + alp1*(Fc(:,T0(:),2).*(ones(T,1)*H2(:)'))...
    + alp1*(Fc(:,T0(:),3).*(ones(T,1)*H3(:)'));

X1=alp1*B/T;% expected background levels 
X=X0+X1;
Fc = Fc.*alp1;
Y=(poissrnd(X)'); % Data generation with Poisson noise   

% save exampleCirc.mat Y Fc Nrow Ncol A0 T0 Nbspectr;

save('exampleCirc.mat','Y','Fc','Nrow','Ncol','A0','T0','Nbspectr');



