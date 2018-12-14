function swfcnp=dcp_nanpad_swfc(swfc,wl)

% 
% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% Function to pad a swfc tensor with nan values along the time axis
% in order to aid in plotting time-window-centered swfc, and to simplify
% comparison to other window sizes and windowless methods
% 

dimms=ndims(swfc);

if dimms==2
    vox=size(swfc,2);
    npad=nan(floor(wl/2),vox);
    swfcnp=cat(1,npad,swfc);
end

if dimms>2
    [~,ch1,ch2]=size(swfc);
    npad=nan(floor(wl/2),ch1,ch2);
    swfcnp=cat(1,npad,swfc);
end

end