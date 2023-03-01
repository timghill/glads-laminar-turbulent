para_proj = load('para_laminar_turbulent.mat');

cd ~/subglacial-emulator/glads/seasonal_KAN/
para_home = load('para_subglacial_emulator');

cd ~/projects/def-gflowers/tghill/laminar-turbulent/glads/01_kan_l_forcing_shmip_topo/

para_proj.physical
para_home.physical

tt = 200*86400 + 365*100*86400;
source_proj = para_proj.input.source_term_c(tt);
source_home = para_home.input.source_term_c(tt);


max(source_proj)
max(source_home)

sum(source_proj)
sum(source_home)

