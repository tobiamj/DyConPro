function [NCNC,NCNCp,sNCNC,sNCNCp]=dcp_neurocog_correlations(X,subvars)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

% 1.  X is an unfolded group dFC tensor
% 2.  subvars is a vector of subject characteristics with which to compute the
%   time course of correlation for each connection in the dFC tensor
% 3.  NCNC is the NeuroCognitive Network Correlation matrix; it holds a
%   correlation value in a time x connection 2D matrix
% 4.  NCNCp is the same as NCNC but holds p-values instead of correlation
%   coefficients


[tt,rr,~]=size(X);
NCNC=zeros(tt,rr);
NCNCp=zeros(tt,rr);

for loop1=1:tt
    [rcor,pval]=corr(squeeze(X(loop1,:,:))',subvars);
    NCNC(loop1,:)=rcor;
    NCNCp(loop1,:)=pval;
end

end
