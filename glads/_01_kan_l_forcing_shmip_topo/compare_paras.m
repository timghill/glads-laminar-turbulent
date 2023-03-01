tt_melt = (0:0.5:365)*86400 + 100*86400*365;
dt = tt_melt(2) - tt_melt(1);

addpath(genpath('../data/functions/'))
addpath(genpath('../data/kan_l_melt/'))
addpath(genpath('../data/topo_x_squared_para/'))

para_proj = load('para_seasonal.mat');
pp_proj = para_proj.physical
dmesh_proj = para_proj.dmesh
melt_proj = zeros(dmesh_proj.tri.n_nodes, 1);

for ii=1:length(tt_melt)
    t = tt_melt(ii);
    melt_proj = melt_proj + dt*para_proj.input.source_term_c(t);
end


cd ~/subglacial-emulator/glads/seasonal_KAN/
para_home = load('~/subglacial-emulator/glads/seasonal_KAN/para_seasonal.mat');
addpath(genpath('~/subglacial-emulator/glads/data/AWS_GEUS/'))
addpath(genpath('~/subglacial-emulator/glads/data/topo_x_squared_para/'))

pp_home = para_home.physical
dmesh_home = para_home.dmesh
melt_home = zeros(dmesh_home.tri.n_nodes, 1);

for ii=1:length(tt_melt)
    t = tt_melt(ii);
    melt_home = melt_home + dt*para_home.input.source_term_c(t);
end

cd ~/projects/def-gflowers/tghill/laminar-turbulent/glads/_01_kan_l_forcing_shmip_topo/

melt_delta = melt_proj - melt_home;
max(abs(melt_delta))
sum(abs(melt_delta))
