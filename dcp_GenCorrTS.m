function [TSM,rtsm,ptsm]=dcp_GenCorrTS(n,ts,corlvl,cmchk)
% 
% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% 1.  n is length of time series
% 2.  ts is number of time series ot generate
% 3.  corlvl is magnitude of correlation to produce
% 4.  cmtx is a target correlation matrix, entered as a vector in brackets [x y;y x]
%        the diagonal of this matrix should be equal to 1s
% 5.  cmchk = 1 (default) means to compute the correlations of the new time series and print to screen

if ~exist('cmchk','var')
    cmchk=1;
end

M=randn(n,ts);
cmtx=eye(ts);cmtx(cmtx==0)=corlvl;

% multiply the matrix with an upper triangular matrix obtain by the Cholesky decomposition of the desired correlation matrix R:
TSM=M*chol(cmtx);
if cmchk==1
    [rtsm,ptsm]=corrcoef(TSM);
%     plot(TSM','LineWidth',2)
end

end