%% Compute time-averaged melt using SHMIP PDD model

% Melt parameters
diurnal = false;

%% Use an easy mesh for plots
x = 1:100e3;
xy = [x; x]';
addpath('../../moulins-base/data/topo_x_squared_para/');
z = bed_elevation_para(xy, 0) + ice_thickness_para(xy, 0);

tt_year = 365*86400;
t0 = 0.1*tt_year;
t1 = 0.75*tt_year;
dt = 86400;

inst_melt = shmip_melt(z, t0, diurnal);
if max(inst_melt)>0
    error('Melt > 0 at initial time')
end

melt = 0*xy(:, 1);
t = t0;
t_melt = 0;
while t<t1
    inst_melt = shmip_melt(z, t, diurnal);
    if max(inst_melt)>0
        t_melt = t_melt + dt;
    end
    melt = melt + dt*inst_melt;
    t = t + dt;
end


if max(inst_melt)>0
    error('Melt > 0 at final time')
end

% scatter(xy(:, 1)/1e3, melt)

t = tiledlayout(1, 1);
ax1 = axes(t);

mean_melt = melt./t_melt;
plot(xy(:, 1)/1e3, mean_melt*86400, 'LineWidth', 2)
hold on

t_peak = 0.5*tt_year;
t_10 = acos(-15/16)/2/pi*tt_year;

peak_melt = shmip_melt(z, t_peak, diurnal);

repr_melt = shmip_melt(z, t_10, diurnal);

% plot(xy(:, 1)/1e3, peak_melt*86400, 'LineWidth', 2)
% plot(xy(:, 1)/1e3, repr_melt*86400, 'LineWidth', 2)

% legend({'Mean melt', 'Peak melt', 'T_{0m} = 10^\circ C'})
grid on

xlabel('x (km)')
ylabel('Melt (m w.e. day^{-1})')

xlim([0, 100])

ax2 = axes(t);
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
xlim([0, 100])

elevation_ticks = [0, 200, 400, 600, 800, 1200, 1400, 1600];
tlabels = {};
tlocs = [];
for ii=1:length(elevation_ticks)
    [~, xpos] = min(abs(z - elevation_ticks(ii)));
    tlocs(ii) = xy(xpos, 1);
    tlabels{ii} = num2str(elevation_ticks(ii));
end
xticks(tlocs/1e3)
xticklabels(tlabels)
xlabel('Elevation (m)')
yticks([])

print('SHMIP_mean_melt', '-dpng', '-r600')

%% Now compute mean on dmesh

% Load mesh
mesh_struct = load('../../moulins-base/data/mesh_1/mesh.mat');
meshes = mesh_struct.meshes;
dmesh = meshes{4};
xy = dmesh.tri.nodes;

addpath('../../moulins-base/data/topo_x_squared_para/');
z = bed_elevation_para(xy, 0) + ice_thickness_para(xy, 0);

tt_year = 365*86400;
t0 = 0.1*tt_year;
t1 = 0.75*tt_year;
dt = 86400;

inst_melt = shmip_melt(z, t0, diurnal);
if max(inst_melt)>0
    error('Melt > 0 at initial time')
end

melt = 0*xy(:, 1);
t = t0;
t_melt = 0;
while t<t1
    inst_melt = shmip_melt(z, t, diurnal);
    if max(inst_melt)>0
        t_melt = t_melt + dt;
    end
    melt = melt + dt*inst_melt;
    t = t + dt;
end

if max(inst_melt)>0
    error('Melt > 0 at final time')
end


mean_melt = melt./t_melt;
annual_melt = melt./tt_year;
dlmwrite('SHMIP_mean_melt.txt', mean_melt)
dlmwrite('SHMIP_annual_melt.txt', annual_melt)

