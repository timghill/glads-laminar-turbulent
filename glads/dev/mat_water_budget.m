fname = './RUN/output_001_steady.mat';

set_paths
addpath(genpath('../data/'));
addpath(genpath('/home/tghill/SFU-code/glads/GlaDS-matlab/'));

out = load(fname);

out_pp = pp_do_pp(out);

dmesh = out.para.dmesh;
dis_nodes = pp_discharge_at_node(out_pp.Q, dmesh.tri);

[dis_into_edges, nn] = pp_sheet_discharge_into_edges(out_pp.h_sheets, out_pp.phi, dmesh.tri, out_pp.para.physical);
Dis_into_edges = dis_into_edges;%.*dmesh.tri.edge_length;
sheet_outflow = (squeeze(sum(Dis_into_edges(dmesh.tri.bmark_edge>0,1,:),1)))';

% Check conservation
ii_ref = 3995;
Q_node_in = dis_nodes(1, ii_ref, end);
Q_node_out = dis_nodes(2, ii_ref, end);

neigh_edges = dmesh.tri.connect_edge_inv{ii_ref};
m_c = zeros(length(neigh_edges), 1);
for jj=1:length(neigh_edges)
    sheet_flux_left = dis_into_edges(neigh_edges(jj), 1, end);
    sheet_flux_right = dis_into_edges(neigh_edges(jj), 2, end);
    sheet_flux_diff = sheet_flux_right - sheet_flux_left;
    m_c(jj) = sheet_flux_diff;
end

int_m_c = m_c.*dmesh.tri.edge_length(neigh_edges);

dis_nodes = pp_discharge_at_node(out_pp.Q, dmesh.tri);

dis_channel_inflow_into_bnodes = (squeeze(sum(dis_nodes(1,dmesh.tri.bmark>0,:),2)))';
dis_channel_outflow_of_bnodes = (squeeze(sum(dis_nodes(2,dmesh.tri.bmark>0,:),2)))';
channel_outflow = dis_channel_inflow_into_bnodes - dis_channel_outflow_of_bnodes;

volume_flux_out = channel_outflow + sheet_outflow .* (sheet_outflow>0);

[~, tot_source_s, source_c] = pp_get_source(dmesh, out.para.input, out.para.time.out_t);
tot_source_ts = tot_source_s + sum(source_c,1); % this only works when de-scaled as is the case here.

dVolume = tot_source_ts - volume_flux_out;

timesteps = size(out_pp.h_sheets,2);
h_vol = sum(dmesh.tri.area_nodes.*out_pp.h_sheets, 1);
S_vol = sum(dmesh.tri.edge_length.*out_pp.S_channels,1);

tt = out_pp.para.time.out_t;
dh_vol = gradient(h_vol)./gradient(tt);
dS_vol = gradient(S_vol)./gradient(tt);

p_w = out_pp.phi - out_pp.para.scale.phi*out_pp.fields.nodes.phi_m;
h_en = p_w.*out_pp.fields.nodes.e_v .* out_pp.para.scale.e_v/1000/9.8;

h_en_scaled = h_e(out.phis - out.fields.nodes.phi_m, out.fields.nodes.e_v, out.para.physical.flags);
h_en_phys = h_en_scaled.*out_pp.para.scale.h;

h_en_vol = sum(dmesh.tri.area_nodes.*h_en, 1);
dh_en_vol = gradient(h_en_vol)./gradient(tt);

dphids = out_pp.grad_phi_edge;
Xi_channel = abs(out_pp.Q.*dphids);
qnorm = squeeze(sqrt(sum(out_pp.qs.^2, 2)));

qc = zeros(dmesh.tri.n_edges, length(tt));
for ii=1:dmesh.tri.n_edges
    neigh_el = dmesh.tri.connect_edge_el(ii, :);
    neighs = neigh_el(neigh_el>0);
    n_neighs = length(neighs);
    for jj=1:n_neighs
        qc(ii, :) = qc(ii, :) + qnorm(neighs(jj), :)/n_neighs;
    end
end

Xi_sheet = abs(out_pp.para.physical.l_c.*qc.*dphids);
diss_channel = Xi_channel/1e3/334e3;
diss_sheet = Xi_sheet/1e3/334e3;
diss = diss_channel + diss_sheet;

tt_day = tt./86400/365;

figure
hold on
plot(h_vol)
plot(h_en_vol)
legend('Sheet', 'englacial')

figure
hold on
title('Flux')
plot(tt_day, sheet_outflow)
plot(tt_day, channel_outflow)
plot(tt_day, sheet_outflow + channel_outflow)
legend('sheet flow', 'channel flow', 'sum')

figure
hold on
title('Storage')
plot(tt_day, dh_vol)
plot(tt_day, dS_vol)
plot(tt_day, dh_en_vol)
plot(tt_day, dh_vol + dS_vol + dh_en_vol)
legend('dh', 'Ds', 'dhen', 'dstorage')

dstorage = dh_vol + dS_vol + dh_en_vol;
outflow = sheet_outflow + channel_outflow;
tot_source_ts;

figure
hold on
plot(tot_source_ts)
plot(outflow)
plot(dstorage)
plot(outflow - dstorage)
legend('source', 'outflow', 'storage', 'outflow - storage')

figure
plot(tt_day, (outflow + dstorage) - tot_source_ts)
grid on
xlabel('Year')
ylabel('Flux non-conservation (m^3 s^{-1})')

print('mass_residual', '-dpng', '-r600')

[volume_flux_out, volume_flux_in, volume_gain, Dh_vol, DS_vol, Dtot_vol, total_volume_gained, total_input_volume, total_output_volume] = pp_check_mass_conservation(out.para,out.fields,out.phis,out.h_sheets, out.S_channels, true);

