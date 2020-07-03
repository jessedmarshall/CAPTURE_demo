function analysisstruct = compute_tsne_features(MLmatobj,mocapstruct,analysisparams)
%subselect features to make a tsne map
% ---------------------------
% (C) Jesse D Marshall 2020
%     Harvard University 


analysisstruct=struct();
numfiles = size(MLmatobj,2);
%mocapstruct_reduced_agg = cell(1,max(chgpsin));
tsnefeat_name = cell(0,1);

agg_feat = [];
    a6_agg = [];
    a6_agg_nooffset = [];
    
    filelength_total = 0;
    fileuseind =[];
    
    featnames_use = {'dyadic_spectrograms_score_wl_appendages_euc','dyadic_spectrograms_score_wl_appendages',...
        'appendage_pca_score_lengths','appendage_pca_score','appendage_pca_score_euc'};
    featnumber = [15,15,10,10,10];
    
    
    %% collect features
    agg_feat_indiv = [];
    
    a_vec =cell(1,5);
    for kz=1:numel(featnames_use)
        try
            ahere = MLmatobj.(featnames_use{kz});
        catch ME
            ahere = MLmatobj.(featnames_use{kz});
        end
        a_vec{kz} = ahere{8}(:,1:min(size(ahere{8},2),featnumber(kz)));
        
        for rr = 1:size(  a_vec{kz} ,2)
            tsnefeat_name{end+1} = strcat(featnames_use{kz},'__',num2str(rr));
        end
        
        %issue with a hanging index too long/short
        agg_feat_indiv = cat(2,agg_feat_indiv,a_vec{kz});
    end
    
    
    %% do the move frames etc here to match the normal
    frameshere_temp = MLmatobj.frames_appendage_gps;
    framesuse = intersect(frameshere_temp{8},mocapstruct.move_frames);
    framesuse = framesuse(1:end-1);%in the old code but not sure why
    framesuse_sub = framesuse(1:analysisparams.tsnegranularity:numel(framesuse));
    
    a6 = framesuse_sub;
    
    [~,framesuse_aggfeat] = intersect(frameshere_temp{8},mocapstruct.move_frames);
    framesuse_aggfeat = framesuse_aggfeat(1:end-1);
    framesuse_aggfeat = framesuse_aggfeat(1:analysisparams.tsnegranularity:numel(framesuse));
    
    a6_agg = cat(1,a6_agg,bsxfun(@plus,reshape(a6,[],1),filelength_total));
    a6_agg_nooffset = cat(1,a6_agg_nooffset,reshape(a6,[],1));
    filelength_total = filelength_total+size(mocapstruct.aligned_mean_position,1);
    
    fileuseind = cat(1,fileuseind,1*ones(numel(a6),1));
    
    
    %% some weird bug for one animal plus the initial case
    if numel(agg_feat)==0
        agg_feat = cat(1,agg_feat,agg_feat_indiv(framesuse_aggfeat,:));
    else
        if size(agg_feat_indiv,2) ==size(agg_feat,2)
            agg_feat = cat(1,agg_feat,agg_feat_indiv(framesuse_aggfeat,:));
        else
            fprintf('weird error where one appendage file is the wrong size \n')
            size(agg_feat)
            size(agg_feat_indiv)
            agg_feat = cat(1,agg_feat,zeros(numel(framesuse_aggfeat),size(agg_feat,2)));
        end
    end
    
    analysisstruct.tsnefeat_name = tsnefeat_name;
    
    frames_track = reshape(1:size(agg_feat,1),[],1);
    
    jt_features = agg_feat(frames_track,:);
    analysisstruct.frames_with_good_tracking{1} = a6_agg(frames_track);
    analysisstruct.frames_tracking_appendages = (frames_track);
    analysisstruct.subset_of_points_to_plot_tsne_capped{1} = 1:numel(frames_track);
    analysisstruct.subset_of_points_to_plot_tsne_move{1} = 1:numel(frames_track);
    analysisstruct.condition_inds = ones(1,numel(frames_track));
    a6_agg_nooffset = a6_agg_nooffset(frames_track);
    fileuseind_restriced = fileuseind(frames_track);
    %% save raw and various parameters
    jt_features_raw = agg_feat(frames_track,:);
    jt_features_mean =  nanmean(jt_features,1);
    jt_features_std =  nanstd(jt_features,[],1);
    
    jt_features = bsxfun(@rdivide,bsxfun(@minus,jt_features,...
        nanmean(jt_features,1)),nanstd(jt_features,[],1));
    
    analysisstruct.jt_features =jt_features;
    analysisstruct.jt_features_raw =jt_features_raw;
    analysisstruct.jt_features_mean =jt_features_mean;
    analysisstruct.jt_features_std =jt_features_std;
    
    agg_feat = [];
    
    
    %% load in the mocapstruct    
    fprintf('loading mocap  \n');
    analysisstruct.file_sizes{1} = [];
    
    frames_with_good_tracking_sub = a6_agg_nooffset;
    
    markernames = fieldnames(mocapstruct.markers_preproc);
    aligned_markers_temp =  mocapstruct.markers_aligned_preproc;
    preproc_markers_temp =  mocapstruct.markers_preproc;
    
    if isfield(mocapstruct,'analog')
        analog_markers_temp =  mocapstruct.analog;
    end
    
    analysisstruct.file_sizes{1} = cat(1,analysisstruct.file_sizes{1}, size(aligned_markers_temp,1));
       
    fieldscopy = {'markernames','fps','links','markercolor','modular_cluster_properties','bad_frames_agg'};
    for fhere= 1:numel(fieldscopy)
        mocapstruct_reduced_agg.(fieldscopy{fhere}) = mocapstruct.(fieldscopy{fhere});
    end
    
    if isfield(mocapstruct,'mocapfiletimes' )
        mocapstruct_reduced_agg.('mocapfiletimes') = mocapstruct.('mocapfiletimes');
    else
        mocapstruct_reduced_agg.('mocapfiletimes') = [];
    end
    
    
    for rr = 1:numel(markernames)
        mocapstruct_reduced_agg.markers_preproc.(markernames{rr}) = ...
            preproc_markers_temp.(markernames{rr})(frames_with_good_tracking_sub ,:);
        
        mocapstruct_reduced_agg.markers_aligned_preproc.(markernames{rr}) =...
            aligned_markers_temp.(markernames{rr})(frames_with_good_tracking_sub ,:);
    end
    
    analysisstruct.mocapstruct_reduced_agg{1} = mocapstruct_reduced_agg;


analysisstruct.coarse_annotation_mat = [];




