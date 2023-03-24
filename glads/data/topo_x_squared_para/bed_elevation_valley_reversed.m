function bed = bed_elevation_valley(xy, time)
% out = bed_elevation(xy, time)
%
% return bed elevation. Maintain time dependence since the GlaDS code
% assumes this is a function of time

% Bed elevation parameters
trough_dz = 350;
ridge_dz = 350;
yc = 12.5e3;
const_bed = 350;

x = xy(:, 1);
y = xy(:, 2);

bed_ramp = trough_dz*(x - min(x))./(max(x) - min(x));

ridge_ramp = ridge_dz*abs(x - max(x))./(max(x) - min(x));
bed_trough = ridge_ramp.*(y - yc).^2./(max(abs(y - yc))).^2;

bed = bed_ramp + bed_trough;

end
