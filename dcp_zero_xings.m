function [zxing,zx_out]=dcp_zero_xings(x)

% 
% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% Input:
% 1. x is a signal
% 
% Output:
% 1. zxing is number of zero crossings in signal x
% 2. zx_out is the time points that strattle the zero crosssings
% 
% 


xx=x;
xx(x<0)=-1;
xx(x>0)=1;
zxspre=find(xx(2:end)~=xx(1:end-1));
zxspost=zxspre+1;
zx_out=[zxspre' zxspost'];
zxing=length(zxspost);


end