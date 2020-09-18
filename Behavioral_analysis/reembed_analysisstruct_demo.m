function analysisstruct=reembed_analysisstruct_demo(analysisstruct,reembed_filename,clustering_filename)

fprintf('loading files for reembedding \n')

clustering_struct = load(clustering_filename);
reembed_struct = load(reembed_filename);
%run knn search
fprintf('starting nearest neighbor search \n')
[mIdx_reembed,mD] = knnsearch(reembed_struct.jtfeatures_agg,analysisstruct.jtfeatures_agg(:,:),'K',25);

%% other NN is built in, in case you want other options
fprintf('starting reembed loop \n')
zValues_reembed = zeros(size(mIdx_reembed,2),2);
zValues_reembed_closemedian = zeros(size(mIdx_reembed,2),2);

zValues_reembed_firstnn = zeros(size(mIdx_reembed,2),2);
zValues_reembed_firstfivenn = zeros(size(mIdx_reembed,2),2);

zValues_reembed_median = zeros(size(mIdx_reembed,2),2);
distvariance = zeros(size(mIdx_reembed,2),1);
zValuesHere.Y = reembed_struct.zValues_importance;
distthresh = 3;
for kj=1:size(mIdx_reembed,1)
   % take the first nearest neighbor
    mediandist = zValuesHere.Y(mIdx_reembed(kj,1),:);
    distfrommedian = sqrt(sum((zValuesHere.Y(mIdx_reembed(kj,:),:)-mediandist).^2,2));
    %find nearby points
        gooddist = find(distfrommedian<distthresh);
distvariance(kj) = nanmean(distfrommedian(gooddist));
%take median
        zValues_reembed_closemedian(kj,:) = median(zValuesHere.Y(mIdx_reembed(kj,gooddist),:),1);
        zValues_reembed_firstfivenn(kj,:) = nanmedian(zValuesHere.Y(mIdx_reembed(kj,1:5),:),1);
end
analysisstruct.zValues_reembed = zValues_reembed_closemedian;
analysisstruct.zValues = zValues_reembed_closemedian;
analysisstruct.reembed_flag = 1;
analysisstruct.reembed_date = datetime('today');
analysisstruct.reembed_filename = reembed_filename;


%% 

%% also get the clustering params
analysisstruct.clustering_filename = clustering_filename;
fprintf('running clustering \n')
analysisstruct.reembed_clusters = 1;
clusterfields = fieldnames(clustering_struct);
for mm = 1:numel(clusterfields)
    analysisstruct.(clusterfields{mm}) = clustering_struct.(clusterfields{mm});
end
analysisstruct.params.reembed =1;
    analysisstruct.params.reorder=1;
analysisstruct = compute_analysis_clusters_demo(analysisstruct,analysisstruct.params);



end