ids = [1, 2, 3, 4, 5];
for ii=1:length(ids)
    fname = sprintf('RUN/output_%03d_steady.mat', ids(ii))
    convert_issm_outputs(fname);
end
