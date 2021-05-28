function output = DoubleGeneralizedGammaCDF(m1,m2,p,q,a1,a2,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DoubleGeneralizedGammaCDF.m 
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
% This script returns the analytical expression of the Cumulative
% Distribution Function (CDF) of the Double Generalized Gamma distribution
%% Parameters 
% T: Threshold
% alpha: positive parameter related to the effective number of large-scale cells of the scattering process
% beta: natural number and it stands for the amount of fading parameter
% Mc: number of Monte Carlo iterations
% a1, a2, omega1, omega2: distribution parameters of the Double Generalized Gamma model
% p, q: positive parameters that satisfy p/q = a1/a2
% m1, m2: shaping parameters modeling the severity of fading
% X: small scale fluctuations
% Y: large scale fluctuations
% I: optical irradiance (atmospheric turbulences)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega1 = (gamma(m1)/gamma(m1+1/a1))^a1 *m1;
omega2 = (gamma(m2)/gamma(m2+1/a2))^a2 *m2;
C = p^(m2-1/2) * q^(m1-1/2) * (2*pi)^(1 - (p+q)/2) / (gamma(m1)*gamma(m2));
z = (x.^a2/omega2)^p * m1^q*m2^p/(p^p*q^q*omega1^q);
an=[1];
ap=[];
bm=[Delta(q,m1) Delta(p,m2)];
bq=[0];
M=meijerG( an, ap, bm, bq, z );
output = C*M;

end