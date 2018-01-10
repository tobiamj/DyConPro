function OutTens=dcp_groupten2ten(X)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

% X is a 4D group tensor
% [time,region,region,subject]=size(X)
%   dimensions 2 and 3 are symmetric, dimension 4 is subjects, and dimension 1 is time
% OutTens is the unfolded and stacked triu of X for each subject

[tt,rr,cc,ss]=size(X);

for loop1=1:ss
    OutTens(:,:,loop1)=dcp_ten2mat(X(:,:,:,loop1));    
end

end


