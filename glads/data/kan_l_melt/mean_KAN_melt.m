%% Compute mean melt for each node in mesh

mesh_nr = 4;
lr = -0.005;
DDF = 0.01/86400;

meshes = load('../mesh/mesh.mat');
dmesh = meshes.meshes{mesh_nr};

t0 = 0;
t1 = 365*86400;
tt = 0:3600:t1;

melt = zeros(dmesh.tri.n_nodes, length(tt));

KAN_data = importdata('KAN_L_2014_temp_clipped.txt');
tt_days = KAN_data(:, 1);
sl_temp = KAN_data(:, 2);

tinterp_days = tt/86400;

temp_interp = interp1(tt_days, sl_temp, tinterp_days, 'linear', 0);

% Elevations and distributed temp
xy = dmesh.tri.nodes;
z = 390 + 6*(sqrt(dmesh.tri.nodes(:, 1)+5e3) - sqrt(5e3));
temp_dist = z*lr + temp_interp;
KAN_melt = DDF*max(0, temp_dist);

dt = tt(2) - tt(1);

% Only average over the melt season, where any melt happens
ii_pos = find(temp_interp>0);
ii_min = min(ii_pos); ii_max = max(ii_pos);
DT = 120*86400;

total_melt = dt*sum(KAN_melt, 2);
mean_melt = total_melt/DT;

save('KAN_mean_melt.mat', 'mean_melt')
