function melt=KAN_moulin_steady(time, pin, dmesh, ii_moulin, catchmap, ra)
% melt = KAN_moulin_steady(time, pin, dmesh, ii_moulin, catchmap, ra)
%
% Compute steady moulin inputs based on KAN_L AWS temperatures.
%
% Returns constant moulin inputs for each specified moulin,
% using average KAN_L melt

ramp = max(0, min(time/86400/365/25 - 1, 1));

melt_struct = load('KAN_mean_melt.mat');
mean_melt = melt_struct.mean_melt;

% Compute catchment melt and moulin inputs
area = dmesh.tri.area_nodes;

moulins = zeros(dmesh.tri.n_nodes, 1);
moulins(ii_moulin) = 1;

catch_melt = zeros(size(moulins));

for i=1:length(ii_moulin)
    mask = catchmap==(i-1);
    catch_area = area(mask);
    catch_melt(ii_moulin(i)) = sum(catch_area.*mean_melt(mask));
end

melt = catch_melt.*ramp;
melt(melt<0) = 0;

end

