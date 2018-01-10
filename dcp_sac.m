function [sac,sacp]=dcp_sac(X,Y)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018
% 
% Function [sac,sacp]=dcp_sac(X,Y) computes the correlation of time-dependent
% spatial maps and a static spatial map. A single time-dependent matrix is denoted as
% Xt, and the static map is Y.  It could be used, for example, to compute
% the correlation between a functional connectivity map and a structural
% connectivity map, or between a functional connectivity map and the
% time-averaged functional connectivity map.
% 
% Inputs:
% 1. X is a dFC tensor with time x roi x roi matrices
% 2. Y is a static matrix that is roi x roi
% 
% Outputs:
% 1. sac is the time course of correlation coefficients for Xt and Y
% 2. sacp is the time course of p-values for correlation of Xt and Y
% 


[td,~,~]=size(X);
triuy=triu(Y,1);
triuy=reshape(triuy,1,numel(triuy));
sac=zeros(td,1);
sacp=zeros(td,1);
for tws=1:td
   triux=triu(squeeze(X(tws,:,:)),1);
   triux=reshape(triux,1,numel(triux));
   [rr,pp]=corrcoef(triux,triuy);
   sac(tws,1)=rr(1,2);
   sacp(tws,1)=pp(1,2);        
end


end











