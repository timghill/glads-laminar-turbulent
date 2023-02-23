meshes = load('../mesh/mesh.mat');
dmesh = meshes.meshes{4};

addpath(genpath('~/SFU-code/glads/GlaDS-matlab/'))
addpath(genpath('~/glads/GlaDS-matlab/'))

n_moulin = 68;

moulindata = readmatrix(sprintf('../moulins/moulins_KAN_%03d.txt', n_moulin));
catchmap = readmatrix(sprintf('../moulins/catchment_map_KAN_%03d.txt', n_moulin));
ii_moulin = moulindata(:, 1) + 1;

pin = [];
source_term_s = make_anon_fn('@(xy, time) double(0.01/86400/365 + 0*xy(:, 1));');
source_term_c = make_anon_fn('@(time) double(KAN_moulin_steady(time, pin, dmesh, ii_moulin, catchmap));', pin, dmesh, ii_moulin, catchmap);

pin.bed_elevation = @(xy, time) 350;
pin.ice_thickness = @(xy, time) 40 + 6*(sqrt(xy(:, 1) + 5e3) - sqrt(5e3));

addpath(genpath('../melt/'))
source_term_shmip = make_anon_fn('@(time) double(source_moulin_shmip(time, pin, dmesh, ii_moulin, catchmap));', pin, dmesh, ii_moulin, catchmap);
time = 100*365*86400;
moulin_in = source_term_c(time);
moulin_in = moulin_in(moulin_in>0);
moulin_in

moulin_shmip = source_term_shmip(time);
moulin_shmip = moulin_shmip(moulin_shmip>0);
moulin_shmip

scatter(moulin_shmip, moulin_in)
grid on