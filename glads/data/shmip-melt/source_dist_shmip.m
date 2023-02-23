function melt = source_dist_shmip(xy, time, pin)
% melt = source_dist_shmip(xy, time, pin)
%
% SHMIP distributed system source term with moulins
%
% Distributed melt is computed as a constant and
% uniform basal melt rate of 0.25 cm/a

bed = pin.bed_elevation(xy, time);
melt = 0.01/86400/365*ones(size(bed));
end
