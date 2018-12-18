function [mdl]=dcp_pca_dim_reduc(X,comps)

% 
% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% Usage:
%     Reduce dimensionality to a set of orthogonal components
%     
% Inputs:
% 1. X is time x channel x channel 3-way array
% 2. comps is number of comps to retain
% 
% Output:
% 1. mdl is dimensionality-reduced data set
% 
% 

sinput=size(X);
if length(sinput)>=3
    xrs=dcp_ten2mat(X);
    [~,mdl]=pcares(xrs,comps);
    mdl=dcp_mat2tens(mdl);
elseif length(sinput)==2
    [~,mdl]=pcares(X,comps);
end

end