 reembed_filename = 'X:\Jesse\tsneanalysisstructs\reembedstruct_full_noflat.mat';
      clustering_filename = 'X:\Jesse\tsneanalysisstructs\clusteringstruct_full_noflat.mat';     
analysisstruct = reembed_analysisstruct_demo(analysisstruct,reembed_filename,clustering_filename);

figure(3)
plot(analysisstruct.zValues_reembed(:,1),analysisstruct.zValues_reembed(:,2),'ob','MarkerFaceColor','b')

% %load cluster names
%  analysisstruct.clusternames_file = 'E:\Dropbox (Olveczky)\JesseMarshall\Bence\Videolabeling\Allday_finaltsne_8\WLW_labels_coarser_JDM3.xlsx';
% [analysisstruct.clusternames,analysisstruct.clusternames_fine,analysisstruct.clusternames_pose,...
%     analysisstruct.clusternames_dynamics,analysisstruct.clusternames_quality,...
%   analysisstruct.clusternames_sorted,analysisstruct.clusternames_unsorted] = ...
%   load_cluster_names_ontology(analysisstruct,analysisstruct.clusternames_file);
% 
rearing_frames = find(ismember(analysisstruct.annot_reordered{end}, ...
   (find(cellfun(@numel,strfind(analysisstruct.clusternames_sorted,'Rear'))))));
animate_markers_nonaligned_fullmovie_demo(analysisstruct.mocapstruct_reduced_agg{1},rearing_frames);