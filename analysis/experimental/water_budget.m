scenarios = {'steady', 'seasonal'};
scenario = scenarios{2};
switch scenario
    case 'steady'
        fname = '../../glads/00_synth_forcing/RUN/output_003_steady.nc';
    case 'seasonal'
        fname = '../../glads/00_synth_forcing/RUN/output_003_seasonal.nc';
end

addpath(genpath('../../glads/data/'));
addpath(genpath('~/SFU-code/glads/GlaDS-matlab/'));
addpath('../../glads/00_synth_forcing/');
load_glads_paths;


basalmelt = 100e3*25e3*0.01/86400/365;

phi = ncread(fname, 'phi');
q = ncread(fname, 'qs');

Q = ncread(fname, 'Q');
S = ncread(fname, 'S_channel');

h = ncread(fname, 'h_sheet');

time = ncread(fname, 'time');

% moulin_input = source_term_c(time');
nodes = ncread(fname, 'nodes');


bed = ncread(fname, 'bed');
phi_bed = 1000*9.8*bed;
pw = phi - phi_bed;

meshes = load('../../glads/data/mesh/mesh.mat');
dmesh = meshes.meshes{4};

% MELT FORCING
addpath('../../glads/data/topo_x_squared_para/')
pin.bed_elevation = make_anon_fn('@(xy, time) double(bed_elevation_flat(xy, time))');
pin.ice_thickness = make_anon_fn('@(xy, time) double(ice_thickness_flat(xy, time))');


addpath('../../glads/data/shmip_melt/')
n_moulin = 68;
moulindata = readmatrix(sprintf('../../glads/data/moulins/moulins_%03d.txt', n_moulin));
catchmap = readmatrix(sprintf('../../glads/data/moulins/catchment_map_%03d.txt', n_moulin));
ii_moulin = moulindata(:, 1) + 1;

source_term_s = make_anon_fn('@(xy, time) double(0.01/86400/365 + 0*xy(:, 1));');

% Moulin inputs will need to be adjusted for diurnal simulations
switch scenario
    case 'steady'
        dt = 86400*365*5;
        tindex = 21;
        source_term_c = make_anon_fn('@(time) double(source_moulin_shmip_adj_steady(time, pin, dmesh, ii_moulin, catchmap));', pin, dmesh, ii_moulin, catchmap);
    case 'seasonal'
        dt = 86400;
        tindex = 365 + 190;
        source_term_c = make_anon_fn('@(time) double(source_moulin_shmip_adj_seasonal(time, pin, dmesh, ii_moulin, catchmap));', pin, dmesh, ii_moulin, catchmap);
end


% CONTINUE
moulin_input = [];
for ii=1:length(time)
    tii = time(ii);
    moulin_input = [moulin_input, source_term_c(tii)];
end
tot_moulin_input = sum(moulin_input, 1);

tot_source_input = sum(source_term_s(nodes, 0).*dmesh.tri.area_nodes).*ones(size(time'));

n_edge = dmesh.tri.n_edges;
qout = 0;
flux_elements = [];
for ii=1:n_edge
    bmark = dmesh.tri.bmark_edge(ii);
    edge_len = dmesh.tri.edge_length(ii);
    if bmark==1
        el_index = max(dmesh.tri.connect_edge_el(ii, :));
        flux_elements = [flux_elements, el_index];
        qii = q(el_index, :, :);
%         qnormal = qii*[-1; 0];
        qnormal = -reshape(qii(1, 1, :), [1, size(q, 3)]);
        qout = qout + qnormal*edge_len;
    end
end

htotal = sum(h.*dmesh.tri.area_nodes, 1);
hrate = (htotal(2:end) - htotal(1:end-1))/dt;
hrate = [hrate, 0];

n_node = dmesh.tri.n_nodes;
Qout = 0;
Qarr = sparse(dmesh.tri.n_nodes, dmesh.tri.n_edges);
Qboundary = [];
for ii=1:n_node
    if dmesh.tri.bmark(ii)==1
        edge_indices = dmesh.tri.connect_edge_inv{ii};
        n_edgeii = length(edge_indices);
        for jj=1:n_edgeii
            edge_num = edge_indices(jj);
            Qjj = Q(edge_num, :);
            if dmesh.tri.connect_edge(edge_num, 1)==ii
                edge_sign = -1;
            elseif dmesh.tri.connect_edge(edge_num, 2)==ii
                edge_sign = 1;
            end
            Q_signed = Qjj*edge_sign;
            Qboundary = [Qboundary; Q_signed];
%             Q_signed = Qjj;
            Qout = Qout + Q_signed;
        end
    end
end

kc = ncread(fname, 'para/cond_c');
S_safe = max(S, 1e-16);
phi1 = phi(dmesh.tri.connect_edge(:, 1), :);
phi2 = phi(dmesh.tri.connect_edge(:, 2), :);
dphids_channel = abs(phi2 - phi1)./dmesh.tri.edge_length;
Xi_channel = abs(dphids_channel.*Q);

qs_norm = sqrt(sum(q.^2, 2));
qs_norm = reshape(qs_norm, [size(qs_norm, 1), size(qs_norm, 3)]);

elements = zeros(dmesh.tri.n_elements, 2);
for ii=1:dmesh.tri.n_elements
    nodes = dmesh.tri.nodes(dmesh.tri.connect(ii, :), :);
    elements(ii, :) = mean(nodes, 1);
end


Xi_sheet = zeros(dmesh.tri.n_edges, length(time));
l_c = ncread(fname, 'para/l_c');
k_s = ncread(fname, 'para/cond_s');
h_interp_fun = scatteredInterpolant(dmesh.tri.nodes, h(:, 1), 'nearest');
for ii=1:length(time)
    h_interp_fun.Values = h(:, ii);
    h_edge = h_interp_fun(dmesh.tri.edge_midpoints);
    q_norm_edge = k_s.*h_edge.^(5./4.).*abs(dphids_channel(:, ii)).^(0.5);
%     q_norm_edge = k_s.*h_edge.*abs(dphids_channel(:, ii));
    sheet_melt = abs(l_c*q_norm_edge.*dphids_channel(:, ii));
    Xi_sheet(:, ii) = sheet_melt;
%     ii
end


Xi = Xi_channel + Xi_sheet;

A_moulin = 10;
V_moulin = sum(A_moulin/1000/9.8 * pw(ii_moulin, :), 1);
Q_moulin = (V_moulin(2:end) - V_moulin(1:end-1))/dt;
Q_moulin = [0, Q_moulin];

e_v = 1e-4;
storage = e_v/1000/9.8*pw.*dmesh.tri.area_nodes;
storage = sum(storage, 1);
Qstorage = (storage(2:end) - storage(1:end-1))/dt;
Qstorage = [Qstorage, 0];

Qmelt = sum(Xi.*dmesh.tri.edge_length, 1)/1000/334e3;

exchange = 0;
edge_exchange = zeros(dmesh.tri.n_edges, length(time));
for ii=1:n_edge
    if dmesh.tri.bmark_edge(ii)==0
        node1 = dmesh.tri.connect_edge(ii, 1);
        node2 = dmesh.tri.connect_edge(ii, 2);
    
        el1 = dmesh.tri.connect_edge_el(ii, 1);
        el2 = dmesh.tri.connect_edge_el(ii, 2);

        el1_midpoint = elements(el1, :);
        el2_midpoint = elements(el2, :);
        d_el1 = dmesh.tri.edge_midpoints(ii, :) - el1_midpoint;
        d_el2 = dmesh.tri.edge_midpoints(ii, :) - el2_midpoint;
    
        edge_vec = dmesh.tri.nodes(node2, :) - dmesh.tri.nodes(node1, :);
        n1 = [edge_vec(2), -edge_vec(1)];

        if dot(n1, d_el1) > 0
            n1 = -n1;
        end

        un1 = n1/sqrt(sum(n1.^2));
        un2 = -un1;

        if dot(un2, d_el2) > 0
            disp('warning')
        end
    
        q1 = q(el1, :, :);
        q2 = q(el2, :, :);
    
        m_c1 = q1(1, 1, :).*un1(1) + q1(1, 2, :).*un1(2);
        m_c2 = q2(1, 1, :).*un2(1) + q2(1, 2, :).*un2(2);
        m_c = reshape(m_c1 + m_c2, [1, length(time)]);
        exchange = exchange + m_c.*dmesh.tri.edge_length(ii);
        edge_exchange(ii, :) = m_c.*dmesh.tri.edge_length(ii);
    end
end
exchange = reshape(exchange, [1, length(time)]);


Svolume = sum(S.*dmesh.tri.edge_length);
Srate = (Svolume(2:end) - Svolume(1:end-1))/dt;
Srate = [0, Srate];

f1 = figure;
plot(qout)
hold on
plot(hrate)
plot(Qout)
plot(Srate)
plot(Qstorage)
plot(Qmelt)
grid on

total = qout + hrate + Qout + Srate + Qstorage + Q_moulin;
% plot(total)
legend({'q', 'dh/dt', 'Q', 'dS/dt', 'e_v', 'Melt opening'}, 'Location', 'north')
print(f1, 'mass_residual_internals', '-dpng', '-r600')

f2 = figure;
hold on
plot(total)
plot(tot_moulin_input + tot_source_input + Qmelt)
legend({'Sum of internal components', 'Melt inputs'}, 'Location', 'northoutside')
grid on
print(f2, 'mass_residual', '-dpng', '-r600')

figure
hold on
plot(total - (tot_moulin_input + tot_source_input + Qmelt))
grid on

% Check mass conservation around nodes
check_nodes = [3995,4125,4012];

for ii=1:dmesh.tri.n_nodes
%     ff = figure;
%     hold on
    nodeindex = ii;
    nodex = dmesh.tri.nodes(nodeindex, 1);
    nodey = dmesh.tri.nodes(nodeindex, 2);
    edge_indices = dmesh.tri.connect_edge_inv{nodeindex};
%     fprintf('Num edges: %d\n', length(edge_indices));
    Q_neigh = Q(edge_indices, :);
    Q_signs = zeros(length(edge_indices), 1);
    tot_edge_melt = 0;
    for jj=1:length(edge_indices)
        neigh_nodes = dmesh.tri.connect_edge(edge_indices(jj), :);
        if neigh_nodes(1)==nodeindex
            Q_signs(jj) = -1;
        elseif neigh_nodes(2)==nodeindex
            Q_signs(jj) = 1;
        end
        
        edgexy = dmesh.tri.edge_midpoints(edge_indices(jj), :);
        Q_signed = Q_neigh(jj, tindex).*Q_signs(jj);
        Qarr(ii, edge_indices(jj)) = Q_signed;
    end
end
%         if ismember(ii, check_nodes)
%         plot(dmesh.tri.nodes(nodeindex, 1), dmesh.tri.nodes(nodeindex, 2), 'ko')
%         if Q_signed>0
%             textcolor = 'k';
%         else
%             textcolor = 'r';
%         end
%         edge_melt = Xi(edge_indices(jj), tindex).*dmesh.tri.edge_length(edge_indices(jj))/334e3/1e3;
%         tot_edge_melt = tot_edge_melt + edge_melt;
%         text(edgexy(1), edgexy(2), num2str(Q_signed), 'Color', textcolor);
%         plot([nodex, edgexy(1)], [nodey, edgexy(2)], 'k')
%         end
%     Qnode = sum(Q_neigh.*Q_signs, 1);
%     exchange_node = sum(edge_exchange(edge_indices, tindex));
%     net = Qnode(tindex) - (exchange_node/2 + tot_edge_melt);
%     tot_edge_melt;
%     title(sprintf('Q net: %.3f, Sum m_c: %.3f', Qnode(tindex), sum(edge_exchange(edge_indices, tindex))))
%     print(ff, sprintf('node_conservation_%d', check_nodes(ii)), '-dpng', '-r600');
% end
% end

Q_cmap = cmocean('turbid');
Qmin = 1; Qmax = 100;

figure('Units', 'inches', 'Position', [2, 2, 8, 4])
T = tiledlayout(1, 1);
ax = nexttile;
hold on
Qplot = abs(Q(:, tindex));
Qplot(Qplot<1) = nan;

Q_sum = full(sum(Qarr, 2));
patch('Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes, 'FaceVertexCData', Q_sum, 'FaceColor', 'flat', 'EdgeColor', 'none')
cmocean('balance')
clim([-1, 1])
axis image
cb = colorbar(ax);
cb.Label.String = 'Q_{\rm{net}} (m^3 s^{-1})';
title('MATLAB-GlaDS')

cax = axes(T);
cax.Visible=false;
edge_plot(gca, dmesh, Qplot, Q_cmap, [1, 100], 'vmin', 1)
colormap(cax, Q_cmap);
cb2 = colorbar(cax);
cb2.Layout.Tile = 'North';
cb2.Label.String = 'Q (m^3 s^{-1})';
% cb2.Ticks = [1, 5, 10, 15, 20];
clim(cax, [Qmin, Qmax])

print('MAT_glads_q_net', '-dpng', '-r600')

