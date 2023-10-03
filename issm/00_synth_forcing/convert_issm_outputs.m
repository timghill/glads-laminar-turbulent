function convert_issm_outputs(mat_fname)
% convert_issm_outputs(fname)
%
% Convert ISSM outputs from *.mat to *.nc

set_paths
nc_fname = replace(mat_fname, '.mat', '.nc');

outs = load(mat_fname);
md = outs.md;

% Output fields
nc.time = [md.results.TransientSolution.time];
nc.N = [md.results.TransientSolution.EffectivePressure];
nc.phi = [md.results.TransientSolution.HydraulicPotential];
nc.h_s = [md.results.TransientSolution.HydrologySheetThickness];
nc.S = [md.results.TransientSolution.ChannelArea];
nc.Q = [md.results.TransientSolution.ChannelDischarge];


% Parameters
nc.k_s = md.hydrology.sheet_conductivity(1);
nc.alpha_s = md.hydrology.sheet_alpha;
nc.beta_s = md.hydrology.sheet_beta;
nc.l_bed = md.hydrology.cavity_spacing;
nc.h_bed = md.hydrology.bump_height(1);
nc.omega = md.hydrology.omega;

nc.k_c = md.hydrology.channel_conductivity;
nc.alpha_c = md.hydrology.channel_alpha;
nc.beta_c = md.hydrology.channel_beta;

nc.e_v = md.hydrology.englacial_void_ratio;

% Constants
nc.rho_w = md.materials.rho_freshwater;
nc.rho_i = md.materials.rho_ice;
nc.mu = md.materials.mu_water;
nc.nu = nc.mu/nc.rho_w;
nc.g = md.constants.g;

% Geometry
nc.bed = md.geometry.bed;
nc.thickness = md.geometry.thickness;
nc.phi_bed = nc.rho_w*nc.g*nc.bed;
nc.phi_0 = nc.phi_bed + nc.rho_i*nc.g*nc.thickness;

% Dimensions
n_nodes = size(nc.phi, 1);
n_edges = size(nc.S, 1);
n_time = size(nc.time, 2);

% Create the file
nccreate(nc_fname, 'time',...
    'Dimensions', {'time', n_time});
ncwrite(nc_fname, 'time', nc.time);

nccreate(nc_fname, 'N',...
    'Dimensions', {'nodes', n_nodes, 'time', n_time});
ncwrite(nc_fname, 'N', nc.N);

nccreate(nc_fname, 'phi',...
    'Dimensions', {'nodes', n_nodes, 'time', n_time});
ncwrite(nc_fname, 'phi', nc.phi);

nccreate(nc_fname, 'h_s',...
    'Dimensions', {'nodes', n_nodes, 'time', n_time});
ncwrite(nc_fname, 'h_s', nc.h_s);

nccreate(nc_fname, 'S',...
    'Dimensions', {'edges', n_edges, 'time', n_time});
ncwrite(nc_fname, 'S', nc.S);

nccreate(nc_fname, 'Q',...
    'Dimensions', {'edges', n_edges, 'time', n_time});
ncwrite(nc_fname, 'Q', nc.Q);

nccreate(nc_fname, 'bed',...
    'Dimensions', {'nodes'});
ncwrite(nc_fname, 'bed', nc.bed);

nccreate(nc_fname, 'thickness', 'Dimensions', {'nodes'});
ncwrite(nc_fname, 'thickness', nc.thickness);

nccreate(nc_fname, 'phi_bed', 'Dimensions', {'nodes'});
ncwrite(nc_fname, 'phi_bed', nc.phi_bed);

nccreate(nc_fname, 'phi_0', 'Dimensions', {'nodes'});
ncwrite(nc_fname, 'phi_0', nc.phi_0);

paras = {'k_s', 'alpha_s', 'beta_s', 'l_bed', 'h_bed', 'omega',...
    'k_c', 'alpha_c', 'beta_c', 'e_v', 'rho_w', 'rho_i', 'mu', 'nu'};

for ii=1:length(paras)
    fieldname = paras{ii};
    nccreate(nc_fname, fieldname);
    ncwrite(nc_fname, fieldname, nc.(fieldname))
end