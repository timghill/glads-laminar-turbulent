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
md.hydrology.sheet_alpha = config.alpha_s;
md.hydrology.sheet_beta = config.beta_s;
md.hydrology.cavity_spacing = config.l_bed;
md.hydrology.bump_height = config.h_bed*ones(md.mesh.numberofvertices, 1);
md.hydrology.channel_sheet_width = config.l_c;
md.hydrology.omega = config.omega;
md.hydrology.englacial_void_ratio = config.e_v;
if md.hydrology.omega > 0
    md.hydrology.istransition = 1;
else
    md.hydrology.istransition = 0;
end

md.hydrology.ischannels = 1;
md.hydrology.channel_conductivity = config.k_c;
md.hydrology.channel_alpha = 5./4.;
md.hydrology.channel_beta = 3./2.;

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
md.hydrology.moulin_input = zeros(md.mesh.numberofvertices, 1);

% Outputs
md.hydrology.requested_outputs = {'default', 'HydrologyWaterVx', 'HydrologyWaterVy'};

% Solve
md.transient = deactivateall(md.transient);
md.transient.ishydrology = 1;

md.cluster = generic('np', 1); % CHANGEME

% Timestepping
hour = 3600;
day = 86400;
md.timestepping.time_step = 6*hour/md.constants.yts;
md.settings.output_frequency = 4*5;

% md.timestepping.final_time = day*30/md.constants.yts;
md.timestepping.final_time = 10*day/md.constants.yts;

% Tolerances
md.stressbalance.restol = 1e-3;
md.stressbalance.reltol = nan;
md.stressbalance.abstol = nan;
md.stressbalance.maxiter = 100;

% Final options
md.verbose.solution = 1;
md.miscellaneous.name = 'seasonal';

end
