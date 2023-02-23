% function melt=KAN_moulin(time, pin, dmesh, ii_moulin, catchmap, ra)

moulins = load('../synthetic_moulins_KAN/moulins_KAN_025.txt');
ii_moulin = moulins(:, 1);
catchmap = load('../synthetic_moulins_KAN/catchment_map_KAN_025.txt');
meshes = load('../mesh/mesh.mat');
dmesh = meshes.meshes{4};

pin.bed_elevation = @(xy, tt) 0*xy(:, 1) + 350;
pin.ice_thickness = @(xy, tt) 6*(sqrt(xy(:, 1) + 5000) - sqrt(5000)) + 390 - 350;

tt = 86400*(0:365) + 5*365*86400;

melt = zeros(dmesh.tri.n_nodes, length(tt));
for ii=1:length(tt)
    melt(:, ii) = KAN_moulin(tt(ii), pin, dmesh, ii_moulin, catchmap, 0);
end

t_plt = tt/86400 - tt(1);
m_plt = 1:length(ii_moulin);
[xx, yy] = meshgrid(t_plt, m_plt);

figure
pcolor(xx, yy, melt(ii_moulin, :))
% caxis([0, 10])
shading flat
cb = colorbar;
cb.Label.String = 'Moulin input (m^3 s^{-1})';
xlabel('Day of year')
ylabel('Moulin')
% cmocean('amp')
print('KAN_moulin_inputs', '-dpng', '-r600')

max(melt(ii_moulin, :), [], 2)
