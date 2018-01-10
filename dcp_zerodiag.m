function zdiag=dcp_zerodiag(x)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

[td,rows,cols]=size(x);
if td>1
	for timestep=1:td
		zdiag(timestep,:,:)=x(timestep,:,:)-diag(diag(x(timestep,:,:)));
	end
else
	zdiag=x-diag(diag(x));
end
