function ML_features = compute_appendage_pc_demos(mocapstruct,ML_features,coeffstruct_in,overwrite_coeff)

%% load in the coefficients
%coeffstruct in is the file
if exist(coeffstruct_in,'file') && overwrite_coeff ==0
    try
        coeffstruct = load(coeffstruct_in);
    catch ME
        coeffstruct = load(coeffstruct_in);
    end
else
    coeffstruct = struct();
end


%% specify the specific angles for the different appendages
appendage_names = {'Head','axial','LArm','Rarm','LLeg','RLeg','','globall','trunk'};
ML_features.appendage_names = appendage_names;

appendage_anglegps{8} = fieldnames(ML_features.jointangle_struct);
appendage_segvals{8} = [1:numel(ML_features.all_seglengths)];

ML_features.appendage_anglegps = appendage_anglegps;
ML_features.appendage_segvals =appendage_segvals;
ML_features.appendage_names =appendage_names;

%% angles to include for mai tsne
appendage_gps = 8;
for kk = appendage_gps
    frames_appendage_gps{kk} = mocapstruct.modular_cluster_properties.clipped_index{8}; %hack for now to all be 8 %mocapstruct.move_frames,
end
ML_features.frames_appendage_gps = frames_appendage_gps;
ML_features.appendage_gps =appendage_gps;
%for each group, get PCs of the features, in the correct range of frames.
%Get local pose, local morphology, and pose dynamics
%hipass clip within a range

%% hipass and clip all of the joint angles
params.fps = 300;
params.difforder = 10;
params.medfiltorder = 3;
params.gaussorder = 2.5;

%% get the appendage ja
appendage_anglevals = cell(1,numel(appendage_names));
appendage_explained = cell(1,numel(appendage_names));
appendage_dyadic_spectrograms = cell(1,numel(appendage_names));
COEFFS_appendages = cell(1,numel(appendage_names));

%% also compute the PCS for the lengths
appendage_lengths = cell(1,numel(appendage_names));
appendage_explained_lengths = cell(1,numel(appendage_names));
appendage_dyadic_spectrograms_lengths = cell(1,numel(appendage_names));
COEFFS_appendages_lengths = cell(1,numel(appendage_names));

COEFFS_appendages_euc = cell(1,numel(appendage_names));

appendage_gps = 8;%[3:6];
fprintf('starting PCA over appendages \n')
for kk = appendage_gps
    fprintf('group %f \n',kk);
    %% for obtaining the pose we don't use the smoothed and clipped struct
    fieldnames_here = appendage_anglegps{kk};
    
    appendage_anglevals{kk} = zeros( numel(ML_features.jointangle_struct.(fieldnames_here{1})(frames_appendage_gps{kk})),...
        numel(fieldnames_here));
    %%get the angles and run PCA
    for zz =1:numel(fieldnames_here)
        appendage_anglevals{kk}(:,zz) =   ML_features.jointangle_struct.(fieldnames_here{zz})(frames_appendage_gps{kk});
    end
    
    meanval_ja{kk} = nanmean(appendage_anglevals{kk},1);
    
    %subtract mean
    appendage_anglevals{kk} = bsxfun(@minus,appendage_anglevals{kk}, meanval_ja{kk} );
    
    %% load coeffs
    coeffname_appendage =strcat('COEFFS_appendages',num2str(kk));
    explainedname_appendage =strcat('EXPLAINED_appendages',num2str(kk));
    
    if (~isfield(coeffstruct,coeffname_appendage) || overwrite_coeff)
        [COEFFS_appendages{kk}, ~, ~, ~,appendage_explained{kk}] = pca(squeeze(appendage_anglevals{kk}));
        coeffstruct.(coeffname_appendage) = COEFFS_appendages{kk};
        coeffstruct.(explainedname_appendage) = appendage_explained{kk};
    else
        COEFFS_appendages{kk} = coeffstruct.(coeffname_appendage);
        appendage_explained{kk} = coeffstruct.(explainedname_appendage);
    end
    
    %multiply out
    appendage_dyadic_spectrograms{kk} =  appendage_anglevals{kk}*COEFFS_appendages{kk};
    
    
    
    %% do for segments
    appendage_lengths{kk} = zeros( numel(ML_features.jointangle_struct.(fieldnames_here{1})(frames_appendage_gps{kk})),...
        numel(appendage_segvals{kk}));
    for ll =1:numel(appendage_segvals{kk})
        appendage_lengths{kk}(:,ll) = ML_features.all_seglengths{appendage_segvals{kk}(ll)}(frames_appendage_gps{kk});
    end
    
    %subtract mean
    appendage_lengths{kk} = bsxfun(@minus,appendage_lengths{kk},nanmean(appendage_lengths{kk},1));
    
    %% load coeffs
    coeffname_appendage_lengths =strcat('COEFFS_appendages_lengths',num2str(kk));
    explainedname_appendage_lengths =strcat('EXPLAINED_appendages_lengths',num2str(kk));
    
    if (~isfield(coeffstruct,explainedname_appendage_lengths) || overwrite_coeff)
        [COEFFS_appendages_lengths{kk}, ~, ~, ~,appendage_explained_lengths{kk}] = pca(squeeze(appendage_lengths{kk}));
        coeffstruct.(coeffname_appendage_lengths) = COEFFS_appendages_lengths{kk};
        coeffstruct.(explainedname_appendage_lengths) = appendage_explained_lengths{kk};
    else
        COEFFS_appendages_lengths{kk} = coeffstruct.(coeffname_appendage_lengths);
        appendage_explained_lengths{kk} = coeffstruct.(explainedname_appendage_lengths);
    end
    appendage_dyadic_spectrograms_lengths{kk} =  (appendage_lengths{kk}*COEFFS_appendages_lengths{kk});
    
    
    %% get the PCS of the euclidean distances as well
    appendage_euc_vecs{kk} = [];
    for ll =1:numel(appendage_segvals{kk})
        appendage_euc_vecs{kk} = cat(2,appendage_euc_vecs{kk},ML_features.all_segments{appendage_segvals{kk}(ll)}(frames_appendage_gps{kk}',:));
    end
    
    %subtract mean
    meanval_euc{kk} = nanmean(appendage_euc_vecs{kk},1);
    appendage_euc_vecs{kk} = bsxfun(@minus,appendage_euc_vecs{kk},    meanval_euc{kk});
    
    %load coeffs
    coeffname_appendage_euc =strcat('COEFFS_appendages_euc',num2str(kk));
    explainedname_appendage_euc =strcat('EXPLAINED_appendages_euc',num2str(kk));
    
    if (~isfield(coeffstruct,coeffname_appendage_euc) || overwrite_coeff)
        [COEFFS_appendages_euc{kk}, ~, ~, ~,appendage_explained_euc{kk}] = pca(squeeze(appendage_euc_vecs{kk}));
        coeffstruct.(coeffname_appendage_euc) = COEFFS_appendages_euc{kk};
        coeffstruct.(explainedname_appendage_euc) = appendage_explained_euc{kk};
    else
        COEFFS_appendages_euc{kk} = coeffstruct.(coeffname_appendage_euc);
        appendage_explained_euc{kk} = coeffstruct.(explainedname_appendage_euc);
    end
    appendage_dyadic_spectrograms_euc{kk} =  appendage_euc_vecs{kk}*COEFFS_appendages_euc{kk};
    
    
    ML_features.appendage_coeffs{kk} = COEFFS_appendages{kk} ;
    ML_features.appendage_pca_score{kk} = appendage_dyadic_spectrograms{kk} ;
    ML_features.appendage_pca_explained{kk} = appendage_explained{kk} ;
    ML_features.appendage_coeffs_euc{kk} = COEFFS_appendages_euc{kk} ;
    ML_features.appendage_pca_score_euc{kk} = appendage_dyadic_spectrograms_euc{kk} ;
    ML_features.appendage_pca_explained_euc{kk} = appendage_explained_euc{kk} ;
    ML_features.appendage_coeffs_lengths{kk}  = COEFFS_appendages_lengths{kk} ;
    ML_features.appendage_pca_score_lengths{kk}  = appendage_dyadic_spectrograms_lengths{kk} ;
    ML_features.appendage_pca_explained_lengths{kk}  = appendage_explained_lengths{kk} ;
    
end

%% ------------------------
%save coeffs regardless
fprintf('saving appendage coefficients \n')
try
    save(coeffstruct_in,'-struct','coeffstruct','-v7.3')
catch ME
    save(coeffstruct_in,'-struct','coeffstruct','-v7.3')
end

%% apply the PCA to the smoothed dynamics



%% clip out the bad frames, which smooths the breaks between values. The filter order is 0.3 hz.
% these are used for generating the wavelets, to prevent
for kk = appendage_gps
    fprintf('hipass clipping group %f JOINT ANGLe pcS \n',kk);
    fieldnames_here = appendage_anglegps{kk};
    appendage_anglevals = zeros( numel(ML_features.jointangle_struct.(fieldnames_here{1})),...
        numel(appendage_anglegps{kk}));
    for zz =1:numel(fieldnames_here)
        appendage_anglevals(:,zz) =   ML_features.jointangle_struct.(fieldnames_here{zz});
    end
    appendage_anglevals = bsxfun(@minus,appendage_anglevals,meanval_ja{kk})*ML_features.appendage_coeffs{kk};
    appendage_anglevals(find(isnan(appendage_anglevals))) = 0;
    appendage_anglevals(find(isinf(appendage_anglevals))) = 0;
    
    ML_features.appendage_joint_angles_pcs_hipass{kk}=hipass_clip_cell(appendage_anglevals,frames_appendage_gps{kk},params);
end



%% loop over segment vectors and add to the array
for kk = appendage_gps
    fprintf('hipass clipping group %f EUCLIDEAN \n',kk);
    %% get the PCS of the euclidean distances as well
    appendage_euc_vecs = [];
    for ll =1:numel(appendage_segvals{kk})
        appendage_euc_vecs = cat(2,appendage_euc_vecs,ML_features.all_segments{appendage_segvals{kk}(ll)}(:,:));
    end
    appendage_euc_vecs = (appendage_euc_vecs- meanval_euc{kk})*ML_features.appendage_coeffs_euc{kk};
    ML_features.appendage_pca_score_euc_hipassclip{kk}=hipass_clip_cell(appendage_euc_vecs,frames_appendage_gps{kk},params);
    
end


%% in addition to joint angle PCS, get the euclidean/eigenpose PCs as well


%% smooth the dyadic spectrograms as well -- get the pose over a local region
do_timescales =0;
if do_timescales
    timescales = [10,33,100];
    for zz = 1:numel(timescales)
        params.gaussorder = timescales(zz)./2;
        gfilter = fspecial('gaussian',[timescales(zz)*6 1], params.gaussorder); %get out to 3 sigma
        % ML_features.(strcat('ja_dyadic_spectrograms_',num2str(timescales(zz))))= convn(ML_features.ja_dyadic_spectrograms,gfilter,'same');
        
        appendage_pca_score_smoothed = cell(1,numel(appendage_anglegps));
        appendage_pca_score_smoothed_lengths = cell(1,numel(appendage_anglegps));
        
        for kk = 1:numel(appendage_anglegps)
            appendage_pca_score_smoothed{kk} = convn(ML_features.appendage_pca_score{kk},gfilter,'same');
            appendage_pca_score_smoothed_lengths{kk} = convn(ML_features.appendage_pca_score_lengths{kk},gfilter,'same');
            
        end
        ML_features.(strcat('appendage_pca_score',num2str(timescales(zz)))) = appendage_pca_score_smoothed;
        ML_features.(strcat('appendage_pca_score_lengths',num2str(timescales(zz)))) = appendage_pca_score_smoothed_lengths;
    end
end



end
