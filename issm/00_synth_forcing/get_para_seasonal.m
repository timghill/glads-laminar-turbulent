function md = get_para_steady(config)
% get_para_steady(config)

% Start model instance
set_paths;
md = model;

% Load mesh
dmesh = config.dmesh;
md.mesh = mesh2d();
md.mesh.x = dmesh.tri.nodes(:, 1);
md.mesh.y = dmesh.tri.nodes(:, 2);
md.mesh.elements = dmesh.tri.connect;
md = meshconvert(md, md.mesh.elements, md.mesh.x, md.mesh.y);

% Parameters and set glads hydrology
md = setmask(md, '', '');
md = parameterize(md, config.par);

md.hydrology = hydrologyglads();

md.hydrology.sheet_conductivity = config.k_s*ones(md.mesh.numberofvertices, 1);
% md.hydrology.sheet_alpha = config.alpha_s;
% md.hydrology.sheet_beta = config.beta_s;
md.hydrology.cavity_spacing = config.l_bed;
md.hydrology.bump_height = config.h_bed*ones(md.mesh.numberofvertices, 1);
md.hydrology.channel_sheet_width = config.l_c;
% md.hydrology.omega = config.omega;
md.hydrology.englacial_void_ratio = config.e_v;

md.hydrology.ischannels = 1;
md.hydrology.channel_conductivity = config.k_c;
% md.hydrology.channel_alpha = 5./4.;
% md.hydrology.channel_beta = 3./2.;

% Initial conditions
md.initialization.watercolumn = 0.2*md.hydrology.bump_height.*ones(md.mesh.numberofvertices, 1);
md.initialization.channelarea = 0*ones(md.mesh.numberofedges, 1);

phi_bed = md.constants.g*md.materials.rho_freshwater*md.geometry.base;
p_ice = md.constants.g*md.materials.rho_ice*md.geometry.thickness;
md.initialization.hydraulic_potential = phi_bed + p_ice;

md.initialization.vel = 30*ones(md.mesh.numberofvertices, 1);
md.initialization.vx = -30*ones(md.mesh.numberofvertices, 1);
md.initialization.vy = 0*ones(md.mesh.numberofvertices, 1);

% Boundary conditions
md.hydrology.spcphi = NaN(md.mesh.numberofvertices, 1);
pos = find(md.mesh.vertexonboundary & md.mesh.x==min(md.mesh.x));
md.hydrology.spcphi(pos) = phi_bed(pos);

md.hydrology.neumannflux = zeros(md.mesh.numberofelements, 1);

% Forcing
md.hydrology.melt_flag = 1;
md.basalforcings.groundedice_melting_rate = 0.05*ones(md.mesh.numberofvertices, 1);
md.basalforcings.geothermalflux = 0;

n_moulin = 68;
addpath('/home/tghill/SFU-code/laminar-turbulent/glads/data/shmip_melt/');
addpath('/home/tghill/SFU-code/laminar-turbulent/glads/data/moulins/');
moulindata = readmatrix(sprintf('/home/tghill/SFU-code/laminar-turbulent/glads/data/moulins/moulins_%03d.txt', n_moulin));
catchmap = readmatrix(sprintf('/home/tghill/SFU-code/laminar-turbulent/glads/data/moulins/catchment_map_%03d.txt', n_moulin));
ii_moulin = moulindata(:,1) + 1;

addpath(genpath('../../glads/data/'));
pin.bed_elevation = @(xy, t) bed_elevation_flat(xy, t);
pin.ice_thickness = @(xy, t) ice_thickness_flat(xy, t);

source_term_c = @(time) source_moulin_shmip_adj_seasonal(time, pin, dmesh, ii_moulin, catchmap, 0);
tt_melt = 0:(1/365):1;
md.hydrology.moulin_input = zeros(md.mesh.numberofvertices+1, length(tt_melt));
for ii=1:length(tt_melt)
    ti_seconds = tt_melt(ii)*md.constants.yts;
    md.hydrology.moulin_input(1:md.mesh.numberofvertices, ii) = source_term_c(ti_seconds);
end

md.hydrology.moulin_input(end, :) = tt_melt;

% Solve
md.transient = deactivateall(md.transient);
md.transient.ishydrology = 1;

md.cluster = generic('np', 1); % CHANGEME

% Timestepping
hour = 3600;
day = 86400;
dt_hours = 2;
out_freq = (24*5/dt_hours);
md.timestepping.time_step = dt_hours*hour/md.constants.yts;
md.settings.output_frequency = out_freq;

% md.timestepping.final_time = day*30/md.constants.yts;
md.timestepping.final_time = 1;

% Tolerances
md.stressbalance.restol = 1e-3;
md.stressbalance.reltol = nan;
md.stressbalance.abstol = nan;
md.stressbalance.maxiter = 100;

% Final options
md.verbose.solution = 1;
md.miscellaneous.name = 'seasonal';

end
