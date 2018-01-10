function OutTens=dcp_mat2tens(X)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

% X is a 2D (time x regions) matrix to be converted to an upper triangle tensor
% OutTens is the tensorized matrix

[tt,cc]=size(X);
newdim=ceil(sqrt(cc*2));
OutTens=zeros(tt,newdim,newdim);

for tim=1:tt
    b=triu(ones(newdim),1);
    b(b==1)=X(tim,:);
    OutTens(tim,:,:)=b+b';
end

end