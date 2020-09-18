function [agg_mocap_structs,agg_snippetinds,agg_mocap_preproc,agg_frame_inds,agg_file_inds] = collect_mocap_snippets_allzvals_demo(analysisstruct,condhere,params)

%num_snippets = 60;
if nargin <3
    snippet_size = 300;
    snippet_res = 5;
    snippet_frac = 1;
    files_use=99;
else
    snippet_size = params.snippet_size;
    snippet_res = params.snippet_res;
    snippet_frac = params.snippet_frac;
    files_use = params.files_use;
end

%loop over each cluster, collect inds and files for each animal
cluster_inds = cell(numel(analysisstruct.conditionnames),1);
cluster_files = cell(numel(analysisstruct.conditionnames),1);
cluster_filenumbers = cell(numel(analysisstruct.conditionnames),1);
presence_matrix = zeros(numel(analysisstruct.conditionnames),1);

agg_mocap_structs = cell(numel(analysisstruct.conditionnames),1);
agg_snippetinds = cell(numel(analysisstruct.conditionnames),1);
agg_mocap_preproc = cell(numel(analysisstruct.conditionnames),1);
agg_frame_inds = cell(numel(analysisstruct.conditionnames),1);
agg_file_inds = cell(numel(analysisstruct.conditionnames),1);

for kk = condhere%:numel(analysisstruct.conditionnames)
    relative_inds = find(analysisstruct.condition_inds == kk);
    if numel(relative_inds)
        %       presence_matrix(kk,ll) = 1;
        absolute_inds = analysisstruct.frames_with_good_tracking{kk};%(relative_inds);
        %num_exs = min(num_snippets,numel(absolute_inds));
        absolute_inds_restriced =absolute_inds;% sort(randsample(absolute_inds,num_exs),'ASCEND'); %absolute_inds(1:num_exs);
       if numel(analysisstruct.filesizes)>1
        if iscell(analysisstruct.filesizes{kk})
        filelengths=     cat(1,0,cumsum(analysisstruct.filesizes{kk}{1}));
        else
         filelengths=     cat(1,0,cumsum(analysisstruct.filesizes{kk}));
         test{1} =analysisstruct.filesizes{kk};
         analysisstruct.filesizes{kk} = test; 
        end
        
       else
        if iscell(analysisstruct.filesizes{1})
        filelengths=     cat(1,0,cumsum(analysisstruct.filesizes{1}{1}));
        else
         filelengths=     cat(1,0,cumsum(analysisstruct.filesizes{1}));
         test{1} =analysisstruct.filesizes{1};
         analysisstruct.filesizes{1} = test; 
        end  
         analysisstruct.filesizes{kk} = analysisstruct.filesizes{1};
           end
       
           
        %find the appropriate file
        absolute_files_restriced = (absolute_inds_restriced-filelengths')';
        filenumbers = zeros(1,numel(absolute_inds_restriced));
        for jj =1:size(absolute_files_restriced,2)
            filenumbers(jj) = find(absolute_files_restriced(:,jj)>0,1,'last');
        end
        %get indicies modulo the file
        cluster_filenumbers{kk} = (filenumbers);
        cluster_files{kk} = analysisstruct.mocapnames{kk}(filenumbers);
        cluster_inds{kk} =      absolute_files_restriced(sub2ind(size(absolute_files_restriced),filenumbers,1:numel(filenumbers)));
    end
end


snippet_sum = -snippet_size:snippet_res:snippet_size;
% loop over all
tic
for kk = condhere%:numel(analysisstruct.conditionnames)
    unique_files = unique(cat(2,cluster_filenumbers{kk,:}));
    unique_files = unique_files(1:min(files_use,numel(unique_files)));
    %  unique_filenames = ( analysisstruct.mocapnames{kk}(unique_files));
    for jj = unique_files
        aa = (analysisstruct.mocapnames{kk}{jj});
        aligned_markers =   aa.markers_aligned_preproc;
        markers_preproc =   aa.markers_preproc;

        if ~exist('markernames','var')
            markernames = fieldnames(aligned_markers);
        end
        fprintf(' for condition %f for unique file %f of %f \n',kk,jj,numel(unique_files));
        
        %% median filter for kyle
        %get clusters with the same filenames
        clustershere =  find(     cluster_filenumbers{kk} == jj);
        indshere = bsxfun(@plus,cluster_inds{kk}(clustershere),snippet_sum');
        indshere(indshere<1) = 1;
        indshere(indshere>(analysisstruct.filesizes{kk}{1}(jj)-5)) = analysisstruct.filesizes{kk}{1}(jj)-5;
        indshere = reshape(indshere,[],1);
        %mark the boundaries between individual examples
        if numel( agg_snippetinds{kk})
            uniquenum_here = bsxfun(@plus,1:numel(clustershere),max(agg_snippetinds{kk}));
            agg_snippetinds{kk} = cat(1,agg_snippetinds{kk}, reshape(repmat(uniquenum_here,numel(snippet_sum),1),[],1));
        else
            uniquenum_here = 1:numel(clustershere);
            agg_snippetinds{kk} = reshape(repmat(uniquenum_here,numel(snippet_sum),1),[],1);
        end
        
        % clustershere = find(cellfun(@numel,strfind(cluster_files{kk,ll},cluster_files{kk,ll}{jj})));
        
        if isempty(agg_mocap_structs{kk})
            for nn =1:numel(markernames)
                agg_mocap_structs{kk}.(markernames{nn}) = ...
                    aligned_markers.(markernames{nn})(indshere,:);
                agg_mocap_preproc{kk}.(markernames{nn}) = ...
                    markers_preproc.(markernames{nn})(indshere,:);
                        
            end
              agg_frame_inds{kk} = ...
                    indshere;    
              agg_file_inds{kk} = ...
                    jj*ones(numel(indshere),1);
        else
            for nn =1:numel(markernames)
                agg_mocap_structs{kk}.(markernames{nn}) = cat(1, agg_mocap_structs{kk}.(markernames{nn}),...
                    aligned_markers.(markernames{nn})(indshere,:));
                agg_mocap_preproc{kk}.(markernames{nn}) = cat(1, agg_mocap_preproc{kk}.(markernames{nn}),...
                    markers_preproc.(markernames{nn})(indshere,:));
                     
            end
             agg_frame_inds{kk}=  cat(1,agg_frame_inds{kk},...
                    sum(analysisstruct.filesizes{kk}{1}(1:(jj-1)))+indshere);
                agg_file_inds{kk} = cat(1,agg_file_inds{kk},...
                    jj*ones(numel(indshere),1));
        end
    end
end

end