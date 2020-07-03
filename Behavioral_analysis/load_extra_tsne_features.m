
function jt_features = load_extra_tsne_features(mocapstruct,MLmatobj,analysisparams)

%% define feature subset
tsne_features = {'ja_dyadic_spectrograms','appearance_features_agg_score_whitened','pose_score',...
    'spectrogram_pcs_wl_head_angle','spectrogram_pcs_wl_trunk_angle',...
    'absolute_velocity_trunk_abs_100','absolute_std_velocity_trunk_abs_100',...
    'rel_velocity_hipR_abs_100','rel_velocity_hipL_abs_100','rel_std_velocity_hipR_abs_100','rel_std_velocity_hipL_abs_100',...
    'rel_velocity_head_abs_100','rel_std_velocity_head_abs_100','rel_velocity_trunk_abs_100','rel_std_velocity_trunk_abs_100',...
    'rel_velocity_trunk_z_100',...
    'rel_velocity_hipR_abs_300','rel_velocity_hipL_abs_300','rel_std_velocity_hipR_abs_300','rel_std_velocity_hipL_abs_300',...
    'rel_velocity_head_abs_300','rel_std_velocity_head_abs_300','rel_velocity_trunk_abs_300','rel_std_velocity_trunk_abs_300',...
    'rel_velocity_trunk_z_300','absolute_velocity_trunk_abs_300','absolute_std_velocity_trunk_abs_300','rel_velocity_head_z_100','rel_std_velocity_head_z_100'};
num_feat = [10,6,10,15,15];

%% get tsne features
tsnefeatname = cell(0,1);
for ll = 1:numel(tsne_features)
    if numel(num_feat)>=ll
        for mm = 1:num_feat(ll)
            tsnefeatname{numel(tsnefeatname)+1} = strcat(tsne_features{ll},'_',num2str(mm));
        end
    else
        tsnefeatname{numel(tsnefeatname)+1} = tsne_features{ll};
    end
    
end

%% load features
jt_features= [];

frameshere_temp = MLmatobj.framelist_true;
    framesuse = intersect(frameshere_temp,mocapstruct.move_frames);
    framesuse = framesuse(1:end-1);
    framesuse_sub = framesuse(1:analysisparams.tsnegranularity:numel(framesuse));

for mm = 1:numel(tsne_features)
    if mm>numel(num_feat)
        num_feat(mm) = 1;
    end
    evenlyspacedfeature = MLmatobj.(tsne_features{mm})(:,1:num_feat(mm));
    %% limitation of .mat files for spacing
   jt_features = cat(2, jt_features, evenlyspacedfeature(framesuse_sub,:));
end


 jt_features = bsxfun(@rdivide,bsxfun(@minus,jt_features,...
        nanmean(jt_features,1)),nanstd(jt_features,[],1));

end