function Y = denSampling(X,P)

% Compute samples from P
% 
% INPUT:
% X        : Density support
% P        : Probability density function
%
% OUTPUT:
% Y        : Log likelihood
%
% Author: Q.Legros
% Ref: [Q. Legros, S. McLaughlin, Y. Altmann and S. Meignen, "Stochastic EM
%       algorithm for fast analysis of single waveform multi-spectral Lidar 
%       data," 2020 28th European Signal Processing Conference (EUSIPCO), 
%       Amsterdam, 2021, pp. 2413-2417, doi: 10.23919/Eusipco47968.2020.9287414


[N,D]=size(P);
% P = P./(sum(P,2)*ones(1,D));


cum_P = cumsum(P,2);
z = rand(N,1);
A=(z*ones(1,D)<cum_P);
ind=D-sum(A,2)+1;
% keyboard
Y=X(ind);
Y=Y(:);
