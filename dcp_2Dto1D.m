function Out1D=dcp_2Dto1D(x)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

% x is 2D symmetric matrix
% out1D is vectorized upper triangle of x

Out1D=x(find(~tril(ones(size(x)))));


end