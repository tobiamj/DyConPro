function n = dcp_peakfind(x)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

    n=find(diff(diff(x) > 0) < 0);
    u=find(x(n+1) > x(n));
    n(u)=n(u)+1;

end