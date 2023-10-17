# Simulations and analysis for laminar--turbulent flux parameterization

## 'analysis/'
Top-level analysis

## Simulation directories:

### `issm/`:

Configure jobs for running the base scenario with the ISSM implementation of GlaDS.

`defaults.par`: Set default parameters for ISSM simulations

`00_synth_forcing/`: Run synthetic forcing scenario with ISSM model

### `glads/`:

Configure jobs for running with MATLAB GlaDS.

`data/`: Model input data

`get_para.m`: Global default parameter settings

**Main cases**:

`00_synth_forcing/`

`01_kan_forcing/`

**Additional cases**:

`00a_shmip_forcing`: Melt forcing identical to SHMIP case D3

`00b_synth_basalmelt`: Reduced basal melt (from 0.05 m w.e. a-1 to 0.01 m w.e. a-1)

`00c_synth_marine`: Full-flotation boundary condition

`01a_kan_adj_forcing`: KAN forcing with melt volume equal to SHMIP D3

`01b_kan_forcing_ev`: Englacial storage parameter 1e-5 (default: 1e-4)

`01c_kan_diurnal`: KAN forcing with prescribed diurnal melt-rate variations

`01d_kan_basalmelt`: KAN forcing with reduced basal melt (from 0.05 m w.e. a-1 to 0.01  w.e. a-1)

`03a_kan_forcing_trough`, `03b_kan_forcing_valley`, `03c_kan_forcing_trough2`: KAN forcing with modified bed topography

**Supplementary experiments**:

`S00_mesh_refinement`: Steady-state mesh refinement test

`S01_parameter_sensitivity`: Parameter sensitivity experiment

