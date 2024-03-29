function para = get_para_steady(config)
% para = get_para_steady(config)
%
% Set para for steady state run

%% Get defaults and unwrap
para = get_para_topo(config);
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);

%% Time
% pt.end = 20*pp.year;
pt.end   = 100*pp.year;  % end time
pt.out_t = pt.start:5*pp.year:pt.end;

%% Synthetic bed topo
addpath('../data/topo_x_squared_para/')
pin.bed_elevation = make_anon_fn('@(xy, time) double(bed_elevation_trough2(xy, time))');
pin.ice_thickness = make_anon_fn('@(xy, time) double(ice_thickness_trough2(xy, time, pin))', pin);

%% Source functions
n_moulin = config.n_moulin;
moulindata = readmatrix(sprintf('../data/moulins/moulins_%03d.txt', n_moulin));
catchmap = readmatrix(sprintf('../data/moulins/catchment_map_%03d.txt', n_moulin));
ii_moulin = moulindata(:, 1) + 1;

addpath(genpath('../data/kan_l_melt/'))
pin.source_term_c = make_anon_fn('@(time) double(KAN_moulin_steady(time, pin, dmesh, ii_moulin, catchmap));', pin, dmesh, ii_moulin, catchmap);

%% Nondimensionalize and re-wrap
[psp, pst, psmd, psin, mesh] = scale_para(pp, pt, pmd, pin, dmesh, ps);
para = wrap_para(pm, pn, pin, ps, pt, pst, psp, pp, mesh, dmesh, psin, pmd, psmd, pcm);
