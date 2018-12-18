function [aic]=dcp_linreg_aic(n,k,y,yhat)

%
% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% Computes the Akaike Information Criteria (AIC) for linear regression model
% 
% Inputs:
% 1. n is sample size (or length of time series y)
% 2. k is number of predictors in the model (DoF)
% 3. y is original data
% 4. yhat is fitted data
% 
% Output:
% 1. AIC
%

k=k+1;
SSE=norm(y-yhat,2)^2;
aic=n*log(SSE/n)+(2*k);

end
