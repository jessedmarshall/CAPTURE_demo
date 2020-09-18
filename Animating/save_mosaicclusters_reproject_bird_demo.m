function save_mosaicclusters_reproject_bird_demo(analysisstruct_temp,conditionhere,videodir,startpt,endpt,ratception_struct)
doalign = 0;
figure(388)
if nargin<4
    endpt = analysisstruct_temp.density_objects;
end
size_thresh = 15; % need at least 15 example points
for ll=startpt:analysisstruct_temp.density_objects
 behinds = find(analysisstruct_temp.annot_reordered{conditionhere}==ll);
 behinds = behinds(1:3:end);

 regioninds = find(ismember(analysisstruct_temp.agg_snippetinds{conditionhere},behinds));
 
 
 mocapstruct_temp = analysisstruct_temp.mocapstruct_reduced_agg{conditionhere};
     markernames = mocapstruct_temp.markernames;
for zz = 1:numel(markernames)
mocapstruct_temp.markers_preproc.(markernames{zz}) = analysisstruct_temp.agg_preproc{conditionhere}.(markernames{zz})(regioninds,:);
mocapstruct_temp.markers_aligned_preproc.(markernames{zz}) = analysisstruct_temp.agg_mocap_structs_snippets{conditionhere}.(markernames{zz})(regioninds,:);

end

%% get inds of each snippet
ind_cell = [];
ind_cell_full = [];
sampleset = randsample(1:numel(behinds),min(numel(behinds),60));

    for jj = sampleset
    %    ind_cell{jj} = find(analysisstruct_temp.agg_snippetinds{conditionhere} ==regioninds(jj));
          %      ind_cell_full{jj} = find(ismember(regioninds,find(analysisstruct_temp.agg_snippetinds{conditionhere} ==regioninds(jj))));
      ind_cell{jj} = find(ismember(regioninds,find(analysisstruct_temp.agg_snippetinds{conditionhere} ==behinds(jj))));
                ind_cell_full{jj} = find(ismember(regioninds,find(analysisstruct_temp.agg_snippetinds{conditionhere} ==behinds(jj))));

        ind_cell{jj} =  ind_cell{jj}(1:numel(ind_cell{jj} ));
    end
    
        for jj = sampleset
if numel(ind_cell{jj} ) == 0
    ind_cell{jj} = [];
 ind_cell_full{jj} = [];
end
        end
        
        
        
goodinds_cell = [];
for jj = 1:numel(ind_cell)
    if numel(ind_cell{jj})
        goodinds_cell = cat(1,goodinds_cell,jj);
    end
end
ind_cell = ind_cell(goodinds_cell);

        goodsampleinds = behinds(sampleset);%(goodinds_cell))
% make movies
if numel(ind_cell) >size_thresh
if doalign
pose_mat = [];
for zz = 1:numel(markernames)
    aa = mocapstruct_temp.markers_aligned_preproc.(markernames{zz});
for mmm=1:3
    badinds = find(abs(aa(:,mmm)-nanmean(aa(:,mmm)))>5*nanstd(aa(:,mmm)));
   aa(badinds,mmm)  = nanmean(aa(setxor(1:size(aa,1),badinds),mmm),1);
end
    pose_mat = cat(2,pose_mat,aa);
end

alignedinds = [];
[alignedinds,clusterinds,alignedinds_unaligned ] =align_snippets_2(pose_mat, ind_cell_full,1:5:30);
[~,sortind] = sort(accumarray(clusterinds,1),'DESCEND');

for nz =1:numel(alignedinds_unaligned)
    alignedinds_unaligned{nz}(alignedinds_unaligned{nz}>size(mocapstruct_temp.markers_preproc.HeadF,1))...
        = size(mocapstruct_temp.markers_preproc.HeadF,1);
end

totalinds = [];
if numel(sortind)>=5
for jk = 1:5
indsuse_temp = find(clusterinds == sortind(jk));
totalinds = cat(1,totalinds,randsample(indsuse_temp,min(10,numel(indsuse_temp))));
end
    clf;
    close(388)
    %animate_markers_aligned_fullmovie(mocapstruct_temp,alignedinds{307})
    M=plot_multi_clusters(mocapstruct_temp,alignedinds_unaligned((totalinds(1:min(25,numel(totalinds))))),numel(ind_cell{1}))
    clf;
      v = VideoWriter(strcat(videodir,filesep,'clusteredmovie_',num2str((ll))),'MPEG-4');
    open(v)
    writeVideo(v, M)
    close(v)
end
else
    
    
    
    mocapstruct_temp.predictions = mocapstruct_temp.markers_preproc;
        mocapstruct_temp.cameras = ratception_struct.cameras;
         mocapstruct_temp.fullpath= ratception_struct.fullpath;
         mocapstruct_temp.basefolder= ratception_struct.basefolder;
                  mocapstruct_temp.camuse= ratception_struct.camuse;
goodsampleinds = sort(goodsampleinds);
goodsampleinds = (goodsampleinds(randsample(1:numel(goodsampleinds),min(30,(numel(goodsampleinds))))));

for lk=1:numel(goodsampleinds)
 clustframes=analysisstruct_temp.frames_with_good_tracking{1}(goodsampleinds(lk));
 %(find(analysisstruct.annot_reordered{end}==clustno));
 %% This is in videoframe space
clustframes_good = unique(bsxfun(@plus,floor(clustframes./ratception_struct.sample_factor),...
    (-300)/ratception_struct.sample_factor:floor(20/ratception_struct.sample_factor):(300)/ratception_struct.sample_factor));
clustframes_good(clustframes_good<1) = 1;
clustframes_good(clustframes_good>(floor(size(ratception_struct.aligned_mean_position,1)./ratception_struct.sample_factor)))...
    =(floor(size(ratception_struct.aligned_mean_position,1))./ratception_struct.sample_factor);
 ind_cell_repro{lk} = clustframes_good;
 end
   % mocapstruct_temp.predictions = mocapstruct_temp.markers_preproc;

    %dannce_reprojection(ratception_struct,unique(floor(regioninds(ind_cell{12})./10))')
%    ind_cell_repro = ind_cell;
%     for lk=1:numel(ind_cell)
%         ind_cell_repro{lk} = unique(floor(regioninds(ind_cell{12})./10))';
%     end

    overwrite=0;
    filename = strcat(videodir,filesep,'clusteredmovie_more_',num2str((ll)),'.mp4');
    if ~exist(filename,'file') || overwrite
                  v = VideoWriter(filename,'MPEG-4');
    open(v)
        M=plot_multi_clusters_repro_cont_bird_demo(ratception_struct,ind_cell_repro,numel(ind_cell_repro{1}),v)
  %  writeVideo(v, M)
    close(v)
    end
end
        %nimate_markers_aligned_fullmovie_instances(mocapstruct_temp,ind_cell)
   

end

end