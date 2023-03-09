% Run laminar, turbulent, and transition transient models

%% Laminar models
% run_steady(5e-2, 0, 3, 2, 'laminar_steady.nc');
% run_seasonal(5e-2, 0, 3, 2, 'laminar_steady.nc', 'laminar_seasonal.nc');
%
% Turbulent models
% k_turb = sqrt( 1.97e-6/(1/1000)*5e-2);
k_turb = 5e-3;
run_steady(k_turb, 0, 5/4, 3/2, 'turbulent_winter_sliding.nc');
run_seasonal_winter(k_turb, 0, 5/4, 3/2, 'turbulent_seasonal_sliding.nc');
%

% % Transition models
% run_steady(5e-2, 1/1000, 3, 2, 'transition_steady_1.5.nc');
% run_seasonal(5e-2, 1/1000, 3, 2, 'transition_steady_1.5.nc', 'transition_seasonal_1.5.nc');
