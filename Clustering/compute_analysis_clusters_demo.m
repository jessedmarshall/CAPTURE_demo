function analysisstruct = compute_analysis_clusters_demo(analysisstruct,params)
analysisstruct = cluster_tsne_maps(analysisstruct,analysisstruct.params);
fprintf('finding cluster velocities \n')
analysisstruct = find_cluster_velocities(analysisstruct);


 clusternames = cell(numel(analysisstruct.sorted_clust_ind),1);
for kk = 1:numel(analysisstruct.sorted_clust_ind)
clusternames{kk} = num2str(kk);
end 

fprintf('reordering annotation \n')
    if isfield(params,'reorder') && params.reorder==1
    [analysisstruct] = reorder_annotation(analysisstruct,params);
    end

end