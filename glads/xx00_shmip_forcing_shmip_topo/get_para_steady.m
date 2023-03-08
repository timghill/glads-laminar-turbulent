function para = get_para_steady(config)
% para = get_para_steady(config)
%
% Set para for steady state run

%% Get defaults and unwrap
addpath('../')
para = get_para(config);
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);

%% Time
% pt.end = 20*pp.year;
pt.end   = 100*pp.year;  % end time
pt.out_t = pt.start:5*pp.year:pt.end;

%% Source functions
n_moulin = config.n_moulin;
addpath('../data/shmip_melt/')
moulindata = readmatrix(sprintf('../data/moulins/moulins_normal_%03d.txt', n_moulin));
catchmap = readmatrix(sprintf('../data/moulins/catchment_map_normal_%03d.txt', n_moulin));
ii_moulin = moulindata(:, 1) + 1;

pin.source_term_s = make_anon_fn('@(xy, t) double(0.01/86400/365 + 0*xy(:, 1));');
pin.source_term_c = make_anon_fn('@(t) double(source_moulin_shmip_steady(t, pin, dmesh, ii_moulin, catchmap));', pin, dmesh, ii_moulin, catchmap);

%% Nondimensionalize and re-wrap
[psp, pst, psmd, psin, mesh] = scale_para(pp, pt, pmd, pin, dmesh, ps);
para = wrap_para(pm, pn, pin, ps, pt, pst, psp, pp, mesh, dmesh, psin, pmd, psmd, pcm);