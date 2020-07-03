function analysisstruct_out = find_cluster_velocities(analysisstruct)

velocities_clustertimes = cell(1,(analysisstruct.density_objects));

for ll = 1:numel(analysisstruct.coarse_annotation_mat)
    %% get the velocity sum
    velsum = [];
    for rr = 1:numel(analysisstruct.coarse_annotation_mat{ll})
    velsum = cat(1,velsum,analysisstruct.coarse_annotation_mat{ll}{rr}.veltrace);
        end
if numel(velsum)
    % get the velocity around frames in the cluster
    for nn = 1:(analysisstruct.density_objects)
        %% find the indicies of the ?unsorted clusters
        clusterinds = find(    analysisstruct.annotation_vec{ll,size(analysisstruct.annotation_vec,2)} == nn);
        %% get the average of each frame belonging to the cluster
        if numel(clusterinds)
        adjustedinds = reshape(unique(bsxfun(@plus,analysisstruct.frames_with_good_tracking{ll}(clusterinds),-30:30)),1,[]);
        adjustedinds(adjustedinds>max(analysisstruct.frames_with_good_tracking{ll})) = max(analysisstruct.frames_with_good_tracking{ll});
        adjustedinds(adjustedinds<1) = 1;
        velocities_clustertimes{nn} = cat(1,velocities_clustertimes{nn},...
            nanmean(velsum(adjustedinds)));
        end
    end
end
end
analysisstruct_out = analysisstruct;

for nn = 1:(analysisstruct.density_objects)
    velocities_clustertimes{nn}  = nanmean(   velocities_clustertimes{nn} );
end

analysisstruct_out.velocities_clustertimes = velocities_clustertimes;
end
