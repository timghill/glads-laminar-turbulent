function melt=KAN_dist_steady(time, pin, dmesh)
    % melt = KAN_dist_steady(time, pin, dmesh)
    % computes steady moulin inputs based on KAN_L AWS temperatures.
    %
    % Uses average melt rate from KAN_mountain_glacier_mean_melt.txt
    %
    % See also KAN_mountain_glacier_compute_average_melt

    ramp = max(0, min(time/86400/365/25 - 1, 1));

    steady_melt = readmatrix('KAN_mountain_glacier_mean_melt.txt');

    melt = ramp.*steady_melt;
end

