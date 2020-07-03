function smoothed_struct = smooth_struct(structin,params)

smoothed_struct = structin;
fnames = fieldnames(structin);

for ll = 1:numel(fnames)
    vectorin = real(smoothed_struct.(fnames{ll}));
tracesmooth = medfilt1(vectorin,params.medfiltorder);
gfilter = fspecial('gaussian',[50 1], params.gaussorder);
tracesmoothed = convn(tracesmooth,gfilter,'same');
smoothed_struct.(fnames{ll}) = tracesmoothed;
end
end