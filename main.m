%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main.m 
%
% Created May, 2021
% Elyes Balti
% The University of Texas at Austin
%
% If you use this code or any (modified) part of it in any publication, please cite 
% the following papers: 
%
% E. Balti and B. K. Johnson, 
% "Tractable Approach to MmWaves Cellular Analysis With FSO Backhauling Under Feedback Delay and Hardware Limitations,"
% in IEEE Transactions on Wireless Communications, vol. 19, no. 1, pp. 410-422, Jan. 2020
% 
% E. Balti and M. Guizani, 
% "Mixed RF/FSO Cooperative Relaying Systems With Co-Channel Interference,"
% in IEEE Transactions on Communications, vol. 66, no. 9, pp. 4014-4027, Sept. 2018
%
% Contact email: ebalti@utexas.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description
% This script produces the analytical and Monte Carlo simulations of the Cumulative Distribution Function (CDF)
% of Double Generalized Gamma distribution. This model is used to generate the Free Space Optical (FSO) Irradiance or also termed 
% as the atmospheric turbulences.
%% Parameters 
% T: Threshold
% alpha: positive parameter related to the effective number of large-scale cells of the scattering process
% beta: natural number and it stands for the amount of fading parameter
% Mc: number of Monte Carlo iterations
% a1, a2, omega1, omeg2: distribution parameters of the Double Generalized Gamma model
% p, q: positive parameters that satisfy p/q = a1/a2
% m1, m2: shaping parameters modeling the severity of fading
% X: small scale fluctuations
% Y: large scale fluctuations
% I: optical irradiance (atmospheric turbulences)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc
T = -10:1:10;
T = db2pow(T);
Mc = 1e6;
m1 = 3;
m2 = 2;
p = 2;
q = 1;
a1 = 4;
a2 = 2;

omega1 = (gamma(m1)/gamma(m1+1/a1))^a1 *m1;
omega2 = (gamma(m2)/gamma(m2+1/a2))^a2 *m2;

X = gamrnd(m1,omega1/m1,Mc,1);
Y = gamrnd(m2,omega2/m2,Mc,1);
X = X.^(1/a1);
Y = Y.^(1/a2);
I = X .* Y;


%% Initialize the CDF vectors
CDFa = zeros(length(T),1);
CDFm = zeros(length(T),1);

for ii=1:length(T)
%% Monte Carlo    
    for kk=1:Mc
       if I(kk) < T(ii)
           CDFm(ii) = CDFm(ii) + 1;
       end      
    end
%% Analytical    
CDFa(ii) = DoubleGeneralizedGammaCDF(m1,m2,p,q,a1,a2,T(ii));
end

CDFm = CDFm/Mc;% Averaging over the number of Monte Carlo iterations
T = pow2db(T);

figure; hold on
plot(T,CDFa,'linewidth',2)
plot(T,CDFm,'*','linewidth',2)
xlabel('Threshold (dB)')
ylabel('Cumulative Distribution Function (CDF)')
title('Double Generalized Gamma Distribution')
legend('Analytical','Monte Carlo','location','northwest')
