function melt = source_dist_shmip_steady(time, xy, pin)
    % melt = source_dist_shmip_steady(time, xy, pin)
    % compute steady distributed melt rate using SHMIP melt parameterization
    %
    % Uses 25 year winter steady-state + 25 year linear ramp-up of
    % melt intensity to ensure model stability.

    ramp = max(0, min(time/86400/365/25 - 1, 1));


    % Read steady surface melt
    steady_melt = readmatrix('SHMIP_adj_mean_melt.txt');

    area = dmesh.tri.area_nodes;
    catch_melt = integrate_melt_by_catchment(ii_moulin, catchmap, area, steady_melt);

    melt = catch_melt.*ramp;
end
