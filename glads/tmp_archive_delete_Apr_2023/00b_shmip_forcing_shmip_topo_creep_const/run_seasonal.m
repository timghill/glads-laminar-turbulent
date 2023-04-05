function run_seasonal(ks, omega, alpha, beta, fname_steady, fname_seasonal)
% run_diurnal(gamma, alpha, beta, fname)
%
% Run with diurnal moulin inputs

set_paths;

%% Parameters
mesh_nr = 4;    % High-resolution mesh
n_moulin = 25;
kc = 0.1;
hr = 0.1;

pa = get_para_seasonal(mesh_nr, n_moulin, omega, kc, ks, alpha, beta, hr, fname_steady, fname_seasonal);

%% Run the model
run_model(pa);
