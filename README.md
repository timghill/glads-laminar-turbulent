# GlaDS simulations and analysis for laminar--turbulent flux parameterization

## Structure

### '/'
`analysis/`: top-level analysis

### `glads/`:

`data/`: Model input data
`get_para.m`: Global default parameter settings

Main cases:

`00_shmip_forcing_shmip_topo/`
`01_kan_l_forcing_shmip_topo/`
`02_shmip_forcing_synth_topo/`
`03_kan_l_forcing_synth_topo/`

Supplementary:

`S00_mesh_refinement': SHMIP steady mesh-refinement
`S01_parameter_sensitivity`: Basic one-at-a-time parameter sensitivity experiments, including different scaling for sheet turbulent conductivity.

