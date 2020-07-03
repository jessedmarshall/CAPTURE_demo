function analysisstruct_out = cluster_tsne_maps(analysisstruct,params)
%% code for clustering

analysisstruct_out = analysisstruct;


%% density, watershed, images
tsnehere =analysisstruct.zValues;
density_maps = cell(1,numel(analysisstruct.conditions_to_run)+1);
density_max = params.expansion_factor*max(tsnehere(:));
density_res = params.density_res ;
density_width = params.density_width;
num_conditions_true = numel(unique(analysisstruct.condition_inds));
num_conditions = num_conditions_true+1;

        fprintf('plotting density maps \n');
if ~isfield(params,'reembed') || params.reembed == 0

figure(480)
for ll =1:num_conditions_true
    [xx,yy,density_maps{ll}] = findPointDensity(tsnehere(find(analysisstruct.condition_inds==ll),:),...
        density_width,[density_res density_res],[-density_max density_max]);
    subplot(1,num_conditions_true+1,ll)
    imagesc(flipud(density_maps{ll}))
    %to try [f,xi] = ksdensity(x)
end

[xx,yy,density_jt] = findPointDensity(tsnehere(:,:),density_width,[density_res density_res],[-density_max density_max]);
subplot(1,num_conditions_true+1,num_conditions_true+1)
imagesc(flipud(density_jt))

density_maps{num_conditions_true+1} = density_jt;
analysisstruct_out.density_maps  = density_maps;



%% Get connected components and watershed
density_watersheds = cell(1,num_conditions);
density_cc = cell(1,num_conditions);
density_stats = cell(1,num_conditions);
    density_threshold = params.density_threshold;

for ll = [1:(num_conditions)]
    density_maps{ll}(density_maps{ll}<density_threshold) = 0;
    inv_density_jt = 1./density_maps{ll};
    inv_density_jt(inv_density_jt>(1/density_threshold)) = (1./density_threshold).^2;
    L=watershed(inv_density_jt);
    L(density_maps{ll}==0) = 0;
    density_cc{ll} = bwconncomp(L);
    density_stats{ll} = regionprops(density_cc{ll},'convexhull');
    density_watersheds{ll} = L;
end

figure(482)
imagesc(flipud(density_watersheds{num_conditions}))


% save parameters
analysisstruct_out.density_watersheds = density_watersheds;
analysisstruct_out.density_cc = density_cc;
analysisstruct_out.density_stats = density_stats;

else
    %NB currently reembeding is only setup for one condition
   density_cc =  analysisstruct_out.density_cc;
   density_watersheds = analysisstruct_out.density_watersheds;
density_stats = analysisstruct_out.density_stats;
density_maps  = analysisstruct_out.density_maps;
L = analysisstruct_out.unsorted_watershed;
xx = analysisstruct_out.xx;
yy = analysisstruct_out.yy;
end

%% get points in each cluster
max_objects = max(density_cc{num_conditions}.NumObjects);%,density_cc{2}.NumObjects);

example_inds = cell(1,max_objects);
cluster_num = zeros(1,max_objects);
cluster_mean_length = zeros(4,max_objects);
example_inds_cond = cell(num_conditions,max_objects,num_conditions);
number_cond = zeros(max_objects,num_conditions+1);
example_inds_cond_video = cell(num_conditions,max_objects,num_conditions);
cluster_nums_inclust_indiv  = cell(num_conditions,max_objects,num_conditions);
annotation_vec = cell(num_conditions,num_conditions);
annotation_labels = cell(1,max_objects);
annotation_labels{1} = 'null';

cluster_nums_inclust = cell(1,num_conditions);
cluster_nums_outofclust = cell(1,num_conditions);
cluster_ratios =cell(1,num_conditions);
%% add these to the struct
analysisstruct_out.subset_of_points_to_plot_tsne_capped{num_conditions} = 1:size(tsnehere,1);
analysisstruct_out.frames_with_good_tracking{num_conditions} = 1:size(tsnehere,1);
analysisstruct_out.subset_of_points_to_plot_tsne_move{num_conditions} = 1:size(tsnehere,1);

%% get all points in each region and convert their coordinates to video space if necessary
for cond_select = num_conditions
    %% clustering
    density_objects = density_cc{cond_select}.NumObjects;
    
    L = flipud(density_watersheds{cond_select});
    s = regionprops(L, 'Centroid');
    
    for k = 1:numel(s)
        c = s(k).Centroid;
        text(c(1), c(2), num2str(k), ... %sprintf('%d', integer(ind)),
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle','Color','white');
    end
    hold off
    
    for ll = 1:num_conditions
        annotation_vec{ll,cond_select} = zeros(1,numel(analysisstruct_out.frames_with_good_tracking{ll}));
    end
    
    
    subset_of_points_to_plot =1:size(tsnehere,1);
    condition_labels = cat(1,ones(1,numel(subset_of_points_to_plot)));%,2*ones(1,numel(subset_of_points_to_plot)));
    
    for jj = 1:density_objects
        fprintf('getting inds for density object %f number of clusters %f \n',jj,density_objects);
        [example_inds_x,example_inds_y]  = ind2sub( density_cc{cond_select}.ImageSize,...
            density_cc{cond_select}.PixelIdxList{jj} );
        verts_x = (round(density_stats{cond_select}(jj).ConvexHull(:,1)));
        verts_x(verts_x<1) = 1;
        verts_x(verts_x>numel(xx)) =numel(xx);
        
        verts_y = (round(density_stats{cond_select}(jj).ConvexHull(:,2)));
        verts_y(verts_y<1) = 1;
        verts_y(verts_y>numel(yy)) =numel(yy);
        %% get the points in the watershed polygon
        IN = inpolygon(tsnehere(:,1),tsnehere(:,2),...
            xx(verts_x) ,yy(verts_y));
        
        
        %% assign points to boundaries not initially handled
                example_inds{jj} =  find(IN);
        annotation_labels{jj+1} = strcat('cluster ',num2str(jj));
        
    end
    
    %% assign points not initially put in clusters bc of discretized boundary
    fprintf('reassigning points not put in clusters \n')
    goodpts = cat(1,example_inds{:});
    badpts = setxor(1:size(tsnehere,1),goodpts);
    bad_tsne_x = tsnehere(badpts,1);
     bad_tsne_y = tsnehere(badpts,2);
   
      nnn = density_watersheds{num_conditions};
nnn(nnn>0) = 1;
B = bwboundaries(((nnn)));
overallMinDistance = inf*ones(1,numel(bad_tsne_x));
minddistance_inds = zeros(1,numel(bad_tsne_x));

for kk = 1:numel(B)
     thisBoundaryX = xx(B{kk}(:, 2));
    thisBoundaryY = yy(B{kk}(:,1));
    distances = sqrt((thisBoundaryX - bad_tsne_x).^2 + (thisBoundaryY - bad_tsne_y).^2);
    [minDistance, indexOfMin] = min(distances,[],2);
lessinds = minDistance'<overallMinDistance;
overallMinDistance(find(lessinds)) = minDistance(find(lessinds));
minddistance_inds(lessinds) = kk;
end

for kk = 1:numel(example_inds)
    example_inds{kk} = cat(1,example_inds{kk},badpts(find(minddistance_inds==kk)));
end

        %% get the bout characteristics
    for jj = 1:density_objects
        for ll = 1:num_conditions
            if (ll ~= num_conditions)
                [~,example_inds_cond{cond_select,jj,ll}] = intersect(find(analysisstruct_out.condition_inds==ll),example_inds{jj});
            else
                [~,example_inds_cond{cond_select,jj,ll}] = intersect(find(analysisstruct_out.condition_inds),example_inds{jj});
            end
            annotation_vec{ll,cond_select}((example_inds_cond{cond_select,jj,ll}))=jj;
        end

        if numel(example_inds{jj})
            temp = rectify_inds(unique(bsxfun(@plus,example_inds{jj},-10:10)),max(analysisstruct_out.subset_of_points_to_plot_tsne_capped{cond_select}));
            test = zeros(1,max(analysisstruct_out.subset_of_points_to_plot_tsne_capped{cond_select}));
            test(temp) = 1;
            characteristics = bwconncomp(test);
            cluster_num(jj) = characteristics.NumObjects;
            cluster_mean_length(:,jj) = [characteristics.NumObjects mean(cellfun(@numel,characteristics.PixelIdxList)) ...
                median(cellfun(@numel,characteristics.PixelIdxList)) std(cellfun(@numel,characteristics.PixelIdxList))];
            
            for kk =1:num_conditions
                if numel(example_inds_cond{cond_select,jj,kk})
                    thresh = 2;
                    diff_input = cat(1,0,diff(example_inds_cond{cond_select,jj,kk}));
                    diff_input(diff_input>thresh) = 0;
                    outputstruct = bwconncomp(diff_input);
                    fulloutput = [];
                    for mm = 1:numel(outputstruct.PixelIdxList)
                        fulloutput = cat(2,fulloutput,...
                            analysisstruct_out.subset_of_points_to_plot_tsne_capped{kk}(analysisstruct_out.subset_of_points_to_plot_tsne_capped{kk}(example_inds_cond{cond_select,jj,kk}(outputstruct.PixelIdxList{mm}(1)))):...
                            analysisstruct_out.subset_of_points_to_plot_tsne_capped{kk}(analysisstruct_out.subset_of_points_to_plot_tsne_capped{kk}(example_inds_cond{cond_select,jj,kk}(outputstruct.PixelIdxList{mm}(end)))));
                    end
                    %pick up single instances
                    fulloutput = cat(2,fulloutput,...
                        reshape(analysisstruct_out.subset_of_points_to_plot_tsne_capped{kk}(analysisstruct_out.subset_of_points_to_plot_tsne_capped{kk}(example_inds_cond{cond_select,jj,kk}(diff_input == 0))),1,[]));
                    
                    temp = (unique(bsxfun(@plus,fulloutput',-3:3)));
                    temp(temp<1) = 1;
                    temp(temp>max(analysisstruct_out.subset_of_points_to_plot_tsne_capped{kk})) = max(analysisstruct_out.subset_of_points_to_plot_tsne_capped{kk});
                    example_inds_cond_video{cond_select,jj,kk} = unique(temp);
                    instances = get_indiv_instances( example_inds_cond_video{cond_select,jj,kk} );
                    cluster_nums_inclust_indiv{cond_select,jj,kk} =  numel(instances);
                    cluster_lengths_perbout{cond_select,jj,kk} =  cellfun(@numel,instances);

                else
                    example_inds_cond_video{cond_select,jj,kk} =[];
                end
            end
        end
    end
    % for kk =1:num_conditions
    %aa = numel(get_indiv_instances(   analysisstruct_out.example_inds_cond_video, {14,analysisstruct_out.sorted_clust_ind(62),5});
     %end
    cluster_nums_inclust{cond_select} = cellfun(@numel,example_inds_cond_video(cond_select,:,cond_select));
    %compare to all other clusters
    cluster_nums_outofclust{cond_select} = cellfun(@numel,example_inds_cond_video(cond_select,:,setxor([1:(num_conditions-1)],cond_select)));
    cluster_ratios{cond_select} = cluster_nums_inclust{cond_select}./(cluster_nums_outofclust{cond_select}+cluster_nums_inclust{cond_select});
end
analysisstruct_out.example_inds_cond_video = example_inds_cond_video;
analysisstruct_out.example_inds_cond = example_inds_cond;
analysisstruct_out.cluster_ratios = cluster_ratios;
analysisstruct_out.annotation_vec = annotation_vec;
analysisstruct_out.cluster_nums_inclust_indiv = cluster_nums_inclust_indiv;
analysisstruct_out.cluster_nums_inclust = cluster_nums_inclust;
analysisstruct_out.cluster_lengths_perbout = cluster_lengths_perbout;


%% Initial annotation/transition matrix
%% reorder the clusters
cond_select = num_conditions;
if ~isfield(params,'reembed') || params.reembed == 0

density_objects = density_cc{cond_select}.NumObjects;
var_array = zeros(10,density_objects);
mean_pose = zeros(36,density_objects);
pose_var = zeros(36,density_objects);
frames_agg =cell(1,density_objects);
for ll = 1:density_cc{cond_select}.NumObjects %good_clusters %
    mlist = intersect([1:10,17,18],1:numel(analysisstruct_out.mocapstruct_reduced_agg{1}.markernames));
    allframes = [];
    if (cond_select == num_conditions)
        for kk = 1:num_conditions-1
            allframes_here =  [];
            for jj = mlist
                if numel(example_inds_cond_video{cond_select,ll,kk})
                    allframes_here =  cat(2,allframes_here,...
                        analysisstruct_out.mocapstruct_reduced_agg{kk}.markers_aligned_preproc.(analysisstruct_out.mocapstruct_reduced_agg{kk}.markernames{jj})(example_inds_cond{cond_select,ll,kk},:));
                end
            end
            allframes = cat(1,allframes,allframes_here);
            
        end
    else
        for jj = mlist
            allframes =  cat(2,allframes,...
                analysisstruct_out.mocapstruct_reduced_agg{cond_select}.markers_aligned_preproc.(analysisstruct_out.mocapstruct_reduced_agg{cond_select}.markernames{jj})(example_inds_cond{cond_select,ll,cond_select},:));
        end
    end
    %   frames_agg{ll} = allframes;
    if numel(allframes)
        mean_pose(1:size(allframes,2),ll) = mean(allframes,1);
        pose_var(1:size(allframes,2),ll) = std(allframes,[],1);
    else
        mean_pose(:,ll) = 0;
        pose_var(:,ll) = 0;
    end
    %   [COEFF, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(allframes);
    %  var_array(:,ll) =  (EXPLAINED(1:10))';
end
pose_dist = squareform(pdist(mean_pose','euclidean'));
figure(1112)
Z = linkage(mean_pose','ward','euclidean');
c = cluster(Z,'maxclust',20);
[newind,sorted_clust_ind] = sort(c,'ASCEND');
%crosstab(c,species)
[H,T,outperm] = dendrogram(Z,density_cc{cond_select}.NumObjects);

xlabel('cluster number (original)')


%% if reembeding use the old sorting

analysisstruct_out.dendrogram_linkage = Z;
analysisstruct_out.dendrogram_lines = H;

figure(1111)
imagesc(pose_dist(outperm,outperm))
imagesc(pose_dist(sorted_clust_ind,sorted_clust_ind))

caxis([0 150])
sorted_clust_ind= outperm;
else
    sorted_clust_ind = analysisstruct.sorted_clust_ind;
end

L = zeros(size(density_watersheds{cond_select}));
Lnew = zeros(size(L));
Lnew2 = (zeros(size(L,1),size(L,2),3));
hold on
vv = parula(93);
for k = 1:numel(sorted_clust_ind)
    L(density_cc{cond_select}.PixelIdxList{k}) = k;
end
for k = 1:numel(sorted_clust_ind)
    indstorelabel = find(L ==  sorted_clust_ind(k));%density_cc{cond_select}.PixelIdxList{sorted_clust_ind(k)};
    Lnew(indstorelabel) = (k);
    
end


%% make the visualization comparing the clusters
figure(487)
subplot(1,2,1)
im=imagesc(flipud(Lnew))
colorbar

subplot(1,2,2)
imagesc(flipud(L))
colorbar

ind = 0;
s = regionprops(flipud(L), 'Centroid');
s2 = regionprops(flipud(Lnew), 'Centroid');

for k = 1:numel(sorted_clust_ind')
    
    % get the ind of the sorted cluster
    c = s2((k)).Centroid;
    
    subplot(1,2,1)
    text(c(1), c(2), num2str((k)), ... %sprintf('%d', integer(ind)),
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','Color','white');
    
    c = s((k)).Centroid;
    
    % give it the sorted label
    subplot(1,2,2)
    text(c(1), c(2), num2str((k)), ... %sprintf('%d', integer(ind)),
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','Color','white');
end
hold off
textColor = 'white';

% save parameters of interest
analysisstruct_out.sorted_clust_ind = sorted_clust_ind;
analysisstruct_out.sorted_watershed = Lnew;
analysisstruct_out.unsorted_watershed = L;
analysisstruct_out.density_objects = density_objects;
analysisstruct_out.xx = xx;
analysisstruct_out.yy = yy; 



end