function [ue,le,avg_env]=dcp_get_envelope(X)

% 
% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% Input:
% 1. X is time x channel matrix, or time x channel x channel dFC tensor
%     X can have only one channel if you want
% 
% Output will have same dimensions as input
% 
% If your input is ifc then your output is dcs
% 

if ndims(X)==3
    X=dcp_ten2mat(X);
    nd=1;
else
    nd=0;
end

ue=zeros(size(X));le=zeros(size(X));
for loop1=1:size(X,2)
    ue(:,loop1)=dcp_getspline(X(:,loop1)',[]); 
    le(:,loop1)=-dcp_getspline(-X(:,loop1)',[]); 
end

avg_env=(ue+le)/2;

if nd==1
    avg_env=dcp_mat2tens(avg_env);
    ue=dcp_mat2tens(ue);
    le=dcp_mat2tens(le);
end
        
end
