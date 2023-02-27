function pa = run_job(k_c, bump_scale, alpha, beta, omega, id)

set_paths;
addpath(genpath('../data/functions/'))

fname_steady = sprintf('output_%03d_steady.nc', id);
fname_seasonal = sprintf('output_%03d_seasonal.nc', id);

% Fixed parameters
config.k_s = 0.25;	% Turbulent conductivity
config.l_c = 10;
config.n_moulin = 68;
config.creep_const_soft = 0;
config.mesh_nr = 4;

% Tuning parameters
config.k_c = k_c;
config.h_bed = bump_scale/20;
config.l_bed = bump_scale;
config.alpha = alpha;
config.beta = beta;
config.omega = omega;

if config.alpha<3 && config.omega==0
    % config.k_s = ( (config.h_bed)^(-0.5) * 1/1000/1.79e-6 * (config.k_s^2));
    % omega = 1/1000;
    % k_turb = ( config.k_s * 1.79e-6 / omega / config.h_bed^(3 - 2*config.alpha))^0.5;
    % config.k_s = k_turb
    % config.k_s = 0.006;
    omega = 1/2000;
    h3 = 1.79e-6/(omega)/config.k_s/157;
    k_s = config.k_s * h3^(1 - config.alpha/3) * 157^(2 - 3/2);
    config.k_s = k_s;
end

config

% Case-specific filenames
config.fname_steady = fname_steady;
config.fname_seasonal = fname_seasonal;

% Call GlaDS for each parameter set
para_steady = get_para_steady(config);
para_steady.physical;
pa = para_steady;
% run_model(para_steady);

 
para_seasonal = get_para_seasonal(config);
run_model(para_seasonal);
