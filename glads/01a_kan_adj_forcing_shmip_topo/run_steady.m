function run_steady(ks, omega, alpha, beta, fname)
% run_steady(gamma, alpha, beta, fname)
%
% Run to steady

set_paths;

%% Parameters
mesh_nr = 4;    % High-resolution mesh
n_moulin = 25;
kc = 0.1;
hr = 0.1;

pa = get_para_steady(mesh_nr, n_moulin, omega, kc, ks, alpha, beta, hr, fname);

%% Run the model
run_model(pa);
