function melt=KAN_moulin(time, pin, dmesh, ii_moulin, catchmap, ra)
% melt = source_moulin_shmip(xy, time, pin)
%
% Compute moulin inputs for SHMIP synthetic case based on KAN_L AWS
%
% Returns constant, catchment-dependent moulin inputs for each moulin,
% using SHMIP melt lapse rate with KAN_L temperature timeseries

%% Parameters
lr = -0.005;           % Lapse rate (K/m)
DDF = 0.01/86400;       % Degree day factor (m/K/s)
T_day = 86400;
T_year = 86400*365;
ra = 0;


%% Melt ramp
ramp = max(0, min(time/T_year/25 - 1, 1));
% if time/T_year>=5
%     ramp = 1;
% else
%     ramp = 0;
% end

%% Read AWS data
KAN_data = importdata('KAN_L_2014_temp_clipped.txt');
tt_days = KAN_data(:, 1);
sl_temp = KAN_data(:, 2);

% Interpolate to find instantaneous melt rate
tt_day_resid = mod(time, T_year)/T_day;
temp_interp = interp1(tt_days, sl_temp, tt_day_resid, 'linear', 0);

% Compute elevations and distributed temperature
xy = dmesh.tri.nodes;
z = pin.bed_elevation(xy, 0) + pin.ice_thickness(xy, 0);
temp_dist = z*lr + temp_interp;
KAN_melt = DDF*max(0, temp_dist);

%% Compute catchment melt and moulin inputs
area = dmesh.tri.area_nodes;

moulins = zeros(dmesh.tri.n_nodes, 1);
moulins(ii_moulin) = 1;

catch_melt = zeros(size(moulins));

for i=1:length(ii_moulin)
    mask = catchmap==(i-1);
    catch_area = area(mask);
    catch_melt(ii_moulin(i)) = sum(catch_area.*KAN_melt(mask));
end

diurnal = 1 - ra*sin(2*pi*time/T_day);
melt = catch_melt.*ramp.*diurnal;
melt(melt<0) = 0;

end
