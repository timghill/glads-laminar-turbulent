function para = get_para(config)
% para = get_para_steady(mesh_nr)
%
% Get default parameters for SHMIP ice-sheet margin domain run

addpath(genpath('../data/functions/'))

% Unpack config
mesh_nr = config.mesh_nr;
n_moulin = config.n_moulin;

filename = config.fname_steady;

para = get_default_para();  % Start with model defaults
[pm, pn, pin, ps, pst, psp, mesh, dmesh, pp, pt, psin, pmd, psmd, pcm] = unwrap_all_para(para);

%% Model description
pm.model_run_descript = 'Run to steady';
pm.git_revision_model_runs = strtrim(git('rev-parse --verify HEAD'));
pm.verbosity = 4;   % Lots of output details

%% Model output directories
pm.dir.model_save = ['./', 'RUN', '/'];
pm.save_filename = [pm.dir.model_save, filename];
pm.save_filename_root = '';
pm.save_index_file = 0;

%% Mesh
pm.file.mesh = '../data/mesh/mesh.mat';
mesh_struct = load(pm.file.mesh);
meshes = mesh_struct.meshes;
dmesh = meshes{mesh_nr};

%%  Physical parameters

pp.cond_c = config.k_c;
% pp.flags.max_S = 500; % limit to channel cross-sectional area for stability
pp.flags.max_S = 5e3;

pp.cond_s = config.k_s;
pp.alpha_s = config.alpha;
pp.beta_s =  config.beta;
pp.omega = config.omega;

pp.l_bed = config.l_bed;
pp.l_c = config.l_c;
pp.h_bed = config.h_bed;

pp.creep_const_s_soft = config.creep_const_soft;

e_v = 1e-5;
pin.e_v = make_anon_fn('@(xy) double(0*xy(:,1) + e_v)',e_v);

pp.float_frac = 0; % used below for BC

%% Numerics
pn.zero_channels_on_boundary = 1;
steppers = {'fEuler', 'ode113', 'adaptive', 'ode15i', 'ode15s'};
pn.ts.stepper =  steppers{5};

st = {'ode15s', 'ode23s', 'ode23t', 'odebim'};  % can also use ode23t, ode23s but ode15s is probbaly best
pn.ts.ode15s.stepper = st{1};

%% Domain geometry
addpath('../data/topo_x_squared_para/')
pin.bed_elevation = make_anon_fn('@(xy, time) double(bed_elevation_para(xy, time))');
pin.ice_thickness = make_anon_fn('@(xy, time) double(ice_thickness_para(xy, time))');


%% INPUT FUNCTIONS

% Sliding
%u_bed = config.u_bed/365/86400;
%pin.u_bed = make_anon_fn('@(xy) double(u_bed + 0*xy(:, 1))', u_bed);


% Boundary conditions

% Dirichlet BC for phi: applied at nodes where bmark is odd
pin.bc_dirichlet_phi = make_anon_fn('@(xy, time, bmark, phi_0, phi_m) double(pp.float_frac * (phi_0-phi_m) + phi_m)', pp);

% Flux BC for phi and h_w: i.e. set phi or h_w such that this flux is
% given. Applied at edges where bmark_edge is even
% zero flux:
pin.bc_flux = make_anon_fn('@(xy, time, bmark_edge) double(zeros(sum(~logical(mod(bmark_edge,2)) & bmark_edge>0),1))');

% Initial conditions
hr = pp.h_bed;
pin.ic_h = make_anon_fn('@(xy, time) double(hr/2 - 0.0*xy(:,1)/max(xy(:, 1)))', hr);
pin.ic_S =  make_anon_fn('@(xy, time) double(0 + 0*xy(:,1))');

%% Nondimensionalize and wrap
ps = set_default_scales(ps, pp, dmesh);
[psp, pst, psmd, psin, mesh] = scale_para(pp, pt, pmd, pin, dmesh, ps);

para = wrap_para(pm, pn, pin, ps, pt, pst, psp, pp, mesh, dmesh, psin, pmd, psmd, pcm);