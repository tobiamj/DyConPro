function pca_filt=dcp_pca_filter(x,compsdel,ckeep)

% 
% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% 1. X is the input data matrix; time x channels (voxels)
% 2. compsdel is the number of leading comps to discard
%     for example, if compsdel=2, then pca_filt with be recomposed 
%     from comps 3 thru end.
% 
% Usage:
%   If you want to remove the global mean from your data you could use this
%   program.  Set compsdel to 1 and enter X as your entire fMRI data for a
%   single subejct.
% 


mu=mean(x);
[eigvecs,scores]=pca(x);

Xhat=scores(:,compsdel+1:ckeep)*eigvecs(:,compsdel+1:ckeep)';
pca_filt=bsxfun(@plus,Xhat,mu);


end