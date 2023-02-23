KAN = load('KAN_mean_melt');
KAN_melt = KAN.mean_melt;

addpath('../melt/')
shmip = shmip_melt(0);

meshes = load('../mesh/mesh.mat');
dmesh = meshes.meshes{4};

scatter(dmesh.tri.nodes(:, 1), KAN_melt*86400)
hold on
scatter(dmesh.tri.nodes(:, 1), shmip*86400)

legend({'KAN', 'SHMIP'})
grid on