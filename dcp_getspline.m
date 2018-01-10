function s = dcp_getspline(x,type)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

N = length(x);
p = dcp_peakfind(x);
if type==1;
    s = spline([0 p N+1],[0 x(p) 0],1:N);
else
    s=pchip([0 p N+1],[0 x(p) 0],1:N);
end

end