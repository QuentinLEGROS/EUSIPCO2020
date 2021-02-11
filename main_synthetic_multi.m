clear all                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     % clear all
close all
clc
warning('off','all')

 
%% Initialize data - example
% [Y,Fc,Nrow,Ncol,A0,T0]=initdata();

%% Load data
data=load('exampleCirc.mat'); 
Y=data.Y;
Fc=data.Fc;
Nrow=data.Nrow;
Ncol=data.Ncol;
A0=data.A0;
T0=data.T0;
Nbspectr = size(Fc,3);

%% Initialize Parameters
c=0.05; % Range regularization parameter
step_Nx = 4; % Depth grid subsampling
stepg = 1e-3; % Gradient step for the M-step
si=2; % neighborhood for clustering -> patches of size (2*si+1)x(2*si+1)
nb_label = 4; % number of cluster







%% EM algorithm
tic
[a_out,T_out]=EM_ranging(Y,Fc,Nbspectr,Nrow,Ncol,c,si,step_Nx,stepg,nb_label);
toc
%% Plot reconstructions
% PlotResults(A0,a_out,T0,T_out,Nrow,Ncol)
a_out = a_out';
figure
subplot(3,2,1)
imagesc(A0(:,:,1))
title('A0 1')
subplot(3,2,3)
imagesc(A0(:,:,2))
title('A0 2')
subplot(3,2,5)
imagesc(A0(:,:,3))
title('A0 3')
subplot(3,2,2)
imagesc(reshape(a_out(:,1),Nrow,Ncol))
title('Aout 1')
subplot(3,2,4)
imagesc(reshape(a_out(:,2),Nrow,Ncol))
title('Aout 2')
subplot(3,2,6)
imagesc(reshape(a_out(:,3),Nrow,Ncol))
title('Aout 3')


figure
subplot(1,2,1)
imagesc(T0)
title('GT : T0')
subplot(1,2,2)
imagesc(reshape(T_out,Nrow,Ncol))
title('Tout')

