function bindat=dcp_bnrz(x,threshold,aval)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

if isempty(aval)==1 || ~exist('aval','var')==1
    aval=0;
end

bindat=x;
if aval==0;
    bindat(bindat<threshold)=0;
    bindat(bindat>threshold)=1;
end
if aval==1
    bindat(abs(bindat)<threshold)=0;
    bindat(abs(bindat)>threshold)=1;    
end

end