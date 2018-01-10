function [Out2D]=dcp_1Dto2D(x)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

%  1. x is a vector
%  2. Out2D is symmetric matrix with zeros on the diagonal

% This is especially useful for converting a column vector of spatial
% factor loadings into a matrix that can be plotted or used to make brain
% renderings of networks, etc.

sizediag=ceil(sqrt(2*length(x)));
b=triu(ones(sizediag),1);
b(b==1)=x;
Out2D=b+b';

end