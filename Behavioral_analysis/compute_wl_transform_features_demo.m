function ML_features = compute_wl_transform_features_demo(mocapstruct,ML_features,coeffstruct_in,overwrite_coeff)

if exist(coeffstruct_in,'file')
    try
coeffstruct = load(coeffstruct_in);
    catch ME
        coeffstruct = load(coeffstruct_in);
    end
else
    coeffstruct = struct();
end

%% spectrogram parameters
opts.fps =300./1;
opts.clustering_window = opts.fps./2;
opts.clustering_overlap = opts.fps./4;
opts.numclusters = 100;
opts.lambda = 0.1; % regularization
opts.num = 1; % number modes (spectrograms to find) (usually want full dimension)

%% setup cluster properties
opts.whiten = 0;
opts.frameNormalize = 0;
opts.clustermethod = 'GMM';
% num = pcuse;
opts.ds = 1; % down sampling
opts.samprate = 100;
opts.params = struct;
opts.params.samplingFreq = 100;
opts.params.numPeriods=25; %distinct number of frequencies to use

opts.params.minF = 1; % min freq to analyze
opts.params.maxF = 25; % max freq to analyze
     spacing = 6; %downsample factor
     hipass_val = 0.5;

    
opts.pcuse = 20;
opts.numclusters = 100;
opts.lambda = 0.1; % regularization

        num_spectrogram_pcs= 15;

params.fps = 300;
params.difforder = 10;
params.medfiltorder = 3;
params.gaussorder = 2.5;
%% get the appendage segment length pcs as well
%% also compute the wavelet for the appendages
%restrict to appropriate timerange
COEFFS_feat_wl_appendages = cell(1,numel(ML_features.appendage_anglegps));
dyadic_spectrograms_score_wl_appendages= cell(1,numel(ML_features.appendage_anglegps));
explained_wl_appendages = cell(1,numel(ML_features.appendage_anglegps));

COEFFS_feat_wl_appendages_euc = cell(1,numel(ML_features.appendage_anglegps));
dyadic_spectrograms_score_wl_appendages_euc= cell(1,numel(ML_features.appendage_anglegps));
explained_wl_appendages_euc = cell(1,numel(ML_features.appendage_anglegps));

  dHipass = designfilt('highpassiir', 'FilterOrder', 3, 'HalfPowerFrequency', hipass_val/(params.fps/2), ...
         'DesignMethod', 'butter');
     [f1_hipass,f2_hipass] = tf(dHipass);

         
     
for kk = 8%ML_features.appendage_gps(ML_features.appendage_gps<6)
fprintf('starting group %f joint angles \n',kk)
    agg_features_wl = [];
    if size(ML_features.appendage_joint_angles_pcs_hipass{kk},1)>10*1
    for ll = 1:size(ML_features.appendage_joint_angles_pcs_hipass{kk},2)
        ML_features.appendage_joint_angles_pcs_hipass{kk}(find(isnan(ML_features.appendage_joint_angles_pcs_hipass{kk}(:,ll))),ll) = 0;
        ML_features.appendage_joint_angles_pcs_hipass{kk}(find(isinf(ML_features.appendage_joint_angles_pcs_hipass{kk}(:,ll))),ll) = 0;
        
        frame_fragment = filtfilt(f1_hipass,f2_hipass,ML_features.appendage_joint_angles_pcs_hipass{kk}(:,ll));
        frame_fragment(find(isnan(frame_fragment))) =0;
        frame_fragment(find(isinf(frame_fragment))) =0;
        frame_fragment = real(frame_fragment);
        opts.samprate = opts.fps./spacing; 
        [wavelet_coeffs,w_map,fr_wavelet]= return_wavelets(frame_fragment(1:spacing:end,1),1:size(frame_fragment(1:spacing:end,1),1),opts);
        w_map = w_map+3;
        w_map(w_map<0) = 0;
        agg_features_wl = cat(2,agg_features_wl,w_map);
    end
    
    %% load from file
     %% load coeffs
    coeffname_appendage =strcat('COEFFS_appendages_wl',num2str(kk));
    explainedname_appendage =strcat('EXPLAINED_appendages_wl',num2str(kk));

              if (~isfield(coeffstruct,coeffname_appendage) || overwrite_coeff)
                                    fprintf('OVERWRITING APPENDAGes WL \n')
    [COEFFS_feat_wl_appendages{kk}, ~, ~, ~,explained_wl_appendages{kk}] = pca(squeeze(agg_features_wl));
    coeffstruct.(coeffname_appendage) = COEFFS_feat_wl_appendages{kk};
        coeffstruct.(explainedname_appendage) = explained_wl_appendages{kk};
              else
              COEFFS_feat_wl_appendages{kk} = coeffstruct.(coeffname_appendage);   
              explained_wl_appendages{kk} = coeffstruct.(explainedname_appendage);   
              end
              agg_features_wl = bsxfun(@minus,agg_features_wl,mean(agg_features_wl,1));
                 %multiply out
              dyadic_spectrograms_score_wl_appendages{kk} = agg_features_wl*COEFFS_feat_wl_appendages{kk};
           

    
    %% replicate elements to fill in
    replication_factor_wl = spacing;%ceil(size(frame_fragment,1)./size(dyadic_spectrograms_score_wl_appendages{kk}(:,1),1));
    dyadic_spectrograms_score_wl_appendages{kk} = repelem( dyadic_spectrograms_score_wl_appendages{kk}(:,1:min(num_spectrogram_pcs,...
        size(dyadic_spectrograms_score_wl_appendages{kk},2))),replication_factor_wl,1);
    
    %do a last pruning
    if (size(dyadic_spectrograms_score_wl_appendages{kk},1)<size(frame_fragment,1))
        dyadic_spectrograms_score_wl_appendages{kk} = cat(1,dyadic_spectrograms_score_wl_appendages{kk},zeros(size(frame_fragment,1)...
            -size(dyadic_spectrograms_score_wl_appendages{kk},1),size(dyadic_spectrograms_score_wl_appendages{kk},2)));
    else
        dyadic_spectrograms_score_wl_appendages{kk}((end-(size(dyadic_spectrograms_score_wl_appendages{kk},1)-size(frame_fragment,1) )):end,:) = [];
    end
    
    dyadic_spectrograms_score_wl_appendages{kk} = cat(1,dyadic_spectrograms_score_wl_appendages{kk},zeros(1,size(dyadic_spectrograms_score_wl_appendages{kk},2)));

%% save wavelet coefficients here

%%  ---------------------------------
fprintf('starting group %f euclidean \n',kk)
    agg_features_wl_euc = [];
    for ll = 1:size(ML_features.appendage_pca_score_euc_hipassclip{kk},2)
       ML_features.appendage_pca_score_euc_hipassclip{kk}(find(isnan(ML_features.appendage_pca_score_euc{kk}(:,ll))),ll) = 0;
        ML_features.appendage_pca_score_euc_hipassclip{kk}(find(isinf(ML_features.appendage_pca_score_euc{kk}(:,ll))),ll) = 0;
        
        frame_fragment = filtfilt(f1_hipass,f2_hipass,ML_features.appendage_pca_score_euc_hipassclip{kk}(:,ll));
        frame_fragment(find(isnan(frame_fragment))) =0;
        frame_fragment(find(isinf(frame_fragment))) =0;
        frame_fragment = real(frame_fragment);
        opts.samprate = opts.fps./spacing; 
        [wavelet_coeffs,w_map,fr_wavelet]= return_wavelets(frame_fragment(1:6:end,1),1:size(frame_fragment(1:6:end,1),1),opts);
        w_map = w_map+3;
        w_map(w_map<0) = 0;
        agg_features_wl_euc = cat(2,agg_features_wl_euc,w_map);
    end
    
     coeffname_appendage =strcat('COEFFS_appendages_wl_euc',num2str(kk));
    explainedname_appendage =strcat('EXPLAINED_appendages_wl_euc',num2str(kk));

              if (~isfield(coeffstruct,coeffname_appendage) || overwrite_coeff)
                  fprintf('OVERWRITING APPENDAGes WL \n')
    [COEFFS_feat_wl_appendages_euc{kk}, ~, ~, ~,explained_wl_appendages_euc{kk}] = pca(squeeze(agg_features_wl_euc));
    coeffstruct.(coeffname_appendage) = COEFFS_feat_wl_appendages_euc{kk};
        coeffstruct.(explainedname_appendage) = explained_wl_appendages_euc{kk};
              else
              COEFFS_feat_wl_appendages_euc{kk} = coeffstruct.(coeffname_appendage);   
              explained_wl_appendages_euc{kk} = coeffstruct.(explainedname_appendage);   
              end
              agg_features_wl_euc = bsxfun(@minus,agg_features_wl_euc,mean(agg_features_wl_euc,1));
                 %multiply out
    dyadic_spectrograms_score_wl_appendages_euc{kk} =agg_features_wl_euc*COEFFS_feat_wl_appendages_euc{kk};
    
    % replicate elements to fill in
    replication_factor_wl = spacing;%ceil(size(frame_fragment,1)./size(dyadic_spectrograms_score_wl_appendages{kk}(:,1),1));
    dyadic_spectrograms_score_wl_appendages_euc{kk} = repelem( dyadic_spectrograms_score_wl_appendages_euc{kk}(:,1:...
        min(num_spectrogram_pcs,size(dyadic_spectrograms_score_wl_appendages_euc{kk},2))),replication_factor_wl,1);
    
    %do a last pruning
    if (size(dyadic_spectrograms_score_wl_appendages_euc{kk},1)<size(frame_fragment,1))
        dyadic_spectrograms_score_wl_appendages_euc{kk} = cat(1,dyadic_spectrograms_score_wl_appendages_euc{kk},zeros(size(frame_fragment,1)...
            -size(dyadic_spectrograms_score_wl_appendages_euc{kk},1),size(dyadic_spectrograms_score_wl_appendages_euc{kk},2)));
    else
        dyadic_spectrograms_score_wl_appendages_euc{kk}((end-(size(dyadic_spectrograms_score_wl_appendages_euc{kk},1)-size(frame_fragment,1) )):end,:) = [];
    end
    
    dyadic_spectrograms_score_wl_appendages_euc{kk} = cat(1,dyadic_spectrograms_score_wl_appendages_euc{kk},zeros(1,size(dyadic_spectrograms_score_wl_appendages_euc{kk},2)));

ML_features.COEFFS_feat_wl_appendages_euc{kk} = COEFFS_feat_wl_appendages_euc{kk};
ML_features.dyadic_spectrograms_score_wl_appendages_euc{kk} = dyadic_spectrograms_score_wl_appendages_euc{kk};
ML_features.explained_wl_appendages_euc{kk} = explained_wl_appendages_euc{kk};


ML_features.COEFFS_feat_wl_appendages{kk} = COEFFS_feat_wl_appendages{kk};
ML_features.dyadic_spectrograms_score_wl_appendages{kk} = dyadic_spectrograms_score_wl_appendages{kk};
ML_features.explained_wl_appendages{kk} = explained_wl_appendages{kk};
    else
     
ML_features.COEFFS_feat_wl_appendages_euc{kk} =[];
ML_features.dyadic_spectrograms_score_wl_appendages_euc{kk} =[];
ML_features.explained_wl_appendages_euc{kk} = [];


ML_features.COEFFS_feat_wl_appendages{kk} = [];
ML_features.dyadic_spectrograms_score_wl_appendages{kk} = [];
ML_features.explained_wl_appendages{kk} = [];   
    end
end



%% ------------------------
%save coeffs regardless
fprintf('saving appendage coefficients  !!! WAVELET !!!\n')
try
save(coeffstruct_in,'-struct','coeffstruct','-v7.3')
catch ME
    save(coeffstruct_in,'-struct','coeffstruct','-v7.3')
end

end
