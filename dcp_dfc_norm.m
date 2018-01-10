function [dfc_norm,dfc_norm_mat]=dcp_dfc_norm(dfc)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

% 1. dfc is input dfc time series from a single subject; i.e., a time x roi x roi 3D tensor
% 2. dfc_norm is the output wherefor each time series is normed
% 3. output wil be a 3D tensor; 
% 4. if you want a 2D matrix of the triu then use dfc_norm_mat for your output

for loop1=1:length(dfc)
    x=dcp_ten2mat(dfc);
    dfc_norm_mat(:,loop1)=x(:,loop1)./norm(x(:,loop1));  
end

dfc_norm=dcp_mat2tens(dfc_norm_mat);