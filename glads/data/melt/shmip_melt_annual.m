function melt = shmip_melt(t)
% SHMIP_MELT : load mean melt array
%

fpath = 'SHMIP_annual_melt.txt';

melt = readmatrix(fpath);

end

