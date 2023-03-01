function pa = run_job(k_c, bump_scale, alpha, beta, omega, id)

set_paths;
addpath(genpath('../data/functions/'))

fname_steady = sprintf('output_%03d_steady.nc', id);
fname_seasonal = sprintf('output_%03d_seasonal.nc', id);

% Fixed parameters
config.k_s = 0.1;	% Laminar conductivity
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
config.e_v = 1e-4;

if config.alpha<3 && config.omega==0
    omega = 1/2000;
    h3 = 1.79e-6/omega/config.k_s/157;
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
save('para_steady.mat', '-struct', 'pa');
run_model(para_steady);


para_seasonal = get_para_seasonal(config);
run_model(para_seasonal);

% Plot the melt forcing
pa = para_seasonal;
save('para_seasonal.mat', '-struct', 'pa');
mdata = load('../data/moulins/moulins_068.txt');
ii_moulin = mdata(:, 1);
tt_melt = (0:(1*86400):(365*86400)) + 100*365*86400;
length(tt_melt)
melt = zeros(config.n_moulin, length(tt_melt));
for ii=1:length(tt_melt)
    mii = pa.input.source_term_c(tt_melt(ii));
    melt(:, ii) = mii(ii_moulin+1);
end

dt = tt_melt(2) - tt_melt(1);
area = 100e3 * 25e3;
mean_melt = sum(melt(:))*dt/area

