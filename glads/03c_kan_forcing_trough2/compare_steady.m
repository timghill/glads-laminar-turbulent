addpath('../data/functions/')
addpath('../data/topo_x_squared_para/')

para_proj = load('para_steady.mat');
addpath('../data/kan_l_melt/');
tt = 200*86400 + 365*100*86400;
source_proj = para_proj.input.source_term_c(tt);

cd ~/subglacial-emulator/glads/seasonal_KAN/
addpath('../data/AWS_GEUS/')
addpath('../data/topo_x_squared_para/')
para_home = load('para_steady.mat');


para_proj.physical
para_home.physical

source_home = para_home.input.source_term_c(tt);

cd ~/projects/def-gflowers/tghill/laminar-turbulent/glads/01_kan_l_forcing_shmip_topo/

max(source_proj)
max(source_home)

sum(source_proj)
sum(source_home)


