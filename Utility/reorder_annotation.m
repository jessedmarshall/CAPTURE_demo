function [analysisstruct_out] = reorder_annotation(analysisstruct_out,params)
%% this reorders clusters according to the sorted cluster indicies, sorted by pose similarity instead of position in tsne map 


if nargin<2
   params.reorder=1; 
end
cond_select = size(analysisstruct_out.annotation_vec,2);
numcond =  size(analysisstruct_out.annotation_vec,2)-1;
annot_reordered = cell(1,numcond);
if params.reorder

for kk = 1:cond_select
    [annot_reordered{kk}] = reorder_annotation_vec(analysisstruct_out.annotation_vec{kk,cond_select},analysisstruct_out.sorted_clust_ind);
end

else
    for kk = 1:cond_select
    annot_reordered{kk}=analysisstruct_out.annotation_vec{kk,end};    
    end
    analysisstruct_out.sorted_clust_ind = 1:max(analysisstruct_out.sorted_clust_ind);
end
%% get only unique fields
indtouse = zeros(1,numel(analysisstruct_out.clusternames));
annot_reordered_unique = cell(1,numcond);
for zz = 1:numel(analysisstruct_out.clusternames)
  matchedstrings = find( cellfun(@numel, strfind(analysisstruct_out.clusternames,analysisstruct_out.clusternames{zz})));
    indtouse(zz) = find(matchedstrings,1,'first');
end

[unique_vals,unique_inds] = unique(indtouse);

%% THESE NEED TO BE REORDERED
unique_clusternames = cell(numel(analysisstruct_out.clusternames),1);

for k =1:numel(analysisstruct_out.clusternames)
   unique_clusternames{k} =  analysisstruct_out.clusternames{analysisstruct_out.sorted_clust_ind(k)};

end

for kk = 1:cond_select%analysisstruct_out.conditions_to_run
    annot_reordered_unique{kk} = annot_reordered{kk};
    
    for zz = 1:numel(analysisstruct_out.clusternames)
    end
end

annot_reordered_matched = cell(1,numel(analysisstruct_out.matchedconds));
candidate_features_matched = cell(1,numel(analysisstruct_out.matchedconds));
structoutput = cell(1,numel(analysisstruct_out.matchedconds));
%also get the dotplot features

candidate_features_matched_dotplot = cell(1,numel(analysisstruct_out.matchedconds));
structoutput_dotplot = cell(1,numel(analysisstruct_out.matchedconds));
for mm =1:numel(analysisstruct_out.matchedconds)
    for ll = 1:numel(analysisstruct_out.matchedconds{mm})
annot_reordered_matched{mm} = cat(2,annot_reordered_matched{mm},annot_reordered_unique{analysisstruct_out.matchedconds{mm}(ll)});
if isfield(analysisstruct_out,'candidate_features')
candidate_features_matched{mm} = cat(1,candidate_features_matched{mm},analysisstruct_out.candidate_features{analysisstruct_out.matchedconds{mm}(ll)});
else
   candidate_features_matched{mm} =[];
end
    
if isfield(analysisstruct_out,'candidate_features_dotplot')
candidate_features_matched_dotplot{mm} = cat(1,candidate_features_matched_dotplot{mm},analysisstruct_out.candidate_features_dotplot{analysisstruct_out.matchedconds{mm}(ll)});
end
    end

    if isfield(analysisstruct_out,'candidate_features')
    structoutput{mm} = cell2struct(mat2cell(candidate_features_matched{mm},size(candidate_features_matched{mm},1),ones(1,size(candidate_features_matched{mm},2) )),...
    analysisstruct_out.tsnefeatname ,2);
    else
         structoutput{mm} =[];
    end
      
if isfield(analysisstruct_out,'candidate_features_dotplot')
      if numel(analysisstruct_out.tsnefeat_name_dotplot)==0
          analysisstruct_out.tsnefeat_name_dotplot = analysisstruct_out.tsnefeatname...
              (1:size(analysisstruct_out.jt_features_dotplot,2));
        end
    structoutput_dotplot{mm} = cell2struct(mat2cell(candidate_features_matched_dotplot{mm},...
        size(candidate_features_matched_dotplot{mm},1),...
        ones(1,size(candidate_features_matched_dotplot{mm},2) )),...
    ((analysisstruct_out.tsnefeat_name_dotplot(1:size(candidate_features_matched_dotplot{mm},2)))) ,2);
end
end


%% save to output
if isfield(analysisstruct_out,'candidate_features_dotplot')
analysisstruct_out.structoutput_dotplot = structoutput_dotplot;
analysisstruct_out.candidate_features_matched_dotplot = candidate_features_matched_dotplot;
end

analysisstruct_out.structoutput = structoutput;
analysisstruct_out.annot_reordered_matched = annot_reordered_matched;
analysisstruct_out.candidate_features_matched = candidate_features_matched;
analysisstruct_out.annot_reordered_unique = annot_reordered_unique;
analysisstruct_out.unique_clusternames = unique_clusternames;
analysisstruct_out.annot_reordered = annot_reordered;
analysisstruct_out.unique_clusternames_null =  cat(1,'null',analysisstruct_out.unique_clusternames) ;
end