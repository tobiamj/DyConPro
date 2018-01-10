function OutMat=dcp_ten2mat(X)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

% X is a symmetric tensor of 3 dimensions wherefor 
%   dimensions 2 and 3 are symmetric, and dimension 1 is time
% OutMat is the unfolded triu of X

[tt,rr,~]=size(X);
OutMat=zeros(tt,((rr^2)-rr)/2);

for tim=1:tt
    XX=squeeze(squeeze(X(tim,:,:)));
    OutMat(tim,:)=XX(find(~tril(ones(size(XX)))));
end

end


