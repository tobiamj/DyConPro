function [mirpad_sig,mfront,mrear]=dcp_mirror_pad(x,flag)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01.1 release 6/18/2018

% Inputs:
% 1. x is a signal time series (row vector),
%   or a matrix with time down the rows (set of col vectors)
% 2. flag is a flag (0) for non-mirroring
% 
% NOTES:
%     Operates on single vector or matrix (time x roi)
% 
% Outputs:
% 1. mirpad_sig is the signal x mirror padded (also transposed)
% 

if ~exist('flag','var')
    flag=1;
end

[rows,cols]=size(x);
if cols>=1
    x=x';
end

mirror=rows/2;
mfront=floor(mirror);
mrear=ceil(mirror);

if flag==1
    xfront=fliplr(x(:,1:mfront));
    xrear=fliplr(x(:,mrear+1:end));
end
if flag~=1
    xfront=x(:,1:mfront);
    xrear=x(:,mrear+1:end);    
end

mirpad_sig=[xfront x xrear]';

end


