function run_job(k_c, k_s, alpha_s, beta_s, omega, id)
% run_job(k_c, k_s, alpha_s, beta_s, omega, id)
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

config.e_v = 1e-4;

config.omega = omega;

if config.alpha_s<3 && config.omega==0
    % Compute potential gradient for turbulent conductivity scaling
    p_w_min = 40*910*9.81;
    p_w_max = 1520*910*9.81;
    gradphi = (p_w_max - p_w_min)/100e3;
    omega = 1/2000;
    nu = 1.79e-6;
    h3 = nu/(omega)/config.k_s/gradphi;
    k_s = config.k_s * h3^(1 - config.alpha_s/3) * gradphi^(2 - 3/2);
    config.k_s = k_s;
end


% dmeshfile = '/home/tghill/scratch/laminar-turbulent/glads/data/mesh/mesh.mat';
% dmeshfile = '/home/tghill/SFU-code/laminar-turbulent/glads/data/mesh/mesh.mat';
dmeshfile = '../../glads/data/mesh/mesh.mat';
dmeshes = load(dmeshfile);
dmesh = dmeshes.meshes{4};
config.dmesh = dmesh;

config.par = '../defaults.par';
config.name = sprintf('seasonal_%03d', id);

md = get_para_seasonal(config);
md=solve(md,'Transient');

if ~isfolder('RUN/')
    mkdir('RUN/')
end
config.fname = sprintf('RUN/output_%03d.mat', id);
save(config.fname, 'md');

