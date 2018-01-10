function Out3D=dcp_1Dto3D(x,sizediag)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

%  x is a vector
%  size is symmetric dimension of Out2D
%  Out3D is time series of symmetric matrix with zeros on the diagonal

% DIMS=ndims(x);
% idx=logical(eye(size(x)));
% sizediag=length(find(idx==1));
tsteps=numel(x)/(sizediag^2);

for loop1=1:tsteps
    xx=x(1:sizediag);
    b=triu(ones(sizediag),1);
    b(b==1)=xx;
    x(1:sizediag)=[];
    Out3D(loop1,:,:)=b+b';
end


end