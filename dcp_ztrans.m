function zout=dcp_ztrans(x)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

% Apply Fisher r to z transformation on a tensor, matrix or vector
% zout is z-transformed input

zout=atanh(x);

end