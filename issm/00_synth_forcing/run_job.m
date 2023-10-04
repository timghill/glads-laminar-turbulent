function run_job(k_s, k_c, alpha_s, beta_s, e_v, omega)
% run_job(k_s, k_c, alpha_s, alpha_c, omega)
% Execute steady and seasonal ISSM-GlaDS simulations

config.k_s = k_s;
config.k_c = k_c;

config.alpha_s = alpha_s;
config.beta_s = beta_s;

config.alpha_c = 5./4.;
config.beta_c = 3./2.;

config.l_bed = 10;
config.h_bed = 0.5;
config.l_c = 10;

config.e_v = e_v;

config.omega = omega;

% dmeshfile = '/home/tghill/scratch/laminar-turbulent/glads/data/mesh/mesh.mat';
dmeshfile = '/home/tghill/SFU-code/laminar-turbulent/glads/data/mesh/mesh.mat';
dmeshes = load(dmeshfile);
dmesh = dmeshes.meshes{4};
config.dmesh = dmesh;

config.par = '../defaults.par';

md = get_para_seasonal(config);
md=solve(md,'Transient');

if ~isfolder('RUN/')
    mkdir('RUN/')
end
config.fname = 'RUN/output_001_steady.mat';
save(config.fname, 'md');

% convert_issm_outputs(config.fname);
% delete(config.fname)
