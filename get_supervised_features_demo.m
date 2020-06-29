function MLmatobj = get_supervised_features_demo(mocapstruct,framelist,markersetind,ratname,savename,overwrite_savefile,...
    overwrite_coeff,savefolder,aligncommonrat)


%% local pose/dynamics features
% load saved eigenposture feature coefficients and dynamics coefficients
eigenposture_save_folder = 'Y:\Jesse\Data\Annotation\Eigenposture_coefficients';
condition_mlfeatures_save_folder = 'Y:\Jesse\Data\Annotation\condition_mlfeatures';

if nargin <8
    savefolder = condition_mlfeatures_save_folder;
end

decimation_factor = 3;
%parameters for hipass clip and clustering
params.fps = 300;

%% spectrogram parameters
opts.fps =300./1;
opts.clustering_window = opts.fps./2;
opts.clustering_overlap = opts.fps./4;
num_eigenpcs = 10;%size(COEFF,1);
num_dynamics_pcs = 3;

%% wavelet parameters
% for gmm model parameters
opts.pcuse = 20;
opts.numclusters = 100;
opts.lambda = 0.1; % regularization
if aligncommonrat
    ratname = 'JDM25';
end

%% coefficient files for the PCA
pose_coeff_file = strcat(eigenposture_save_folder,filesep,'pose_coeff_',ratname,'.mat');
dyn_coeff_file = strcat(eigenposture_save_folder,filesep,'dynamics_coeff_',ratname,'.mat');
jadyn_coeff_file = strcat(eigenposture_save_folder,filesep,'ja_coeff_',ratname,'.mat');
ja_pose_coeffs_file = strcat(eigenposture_save_folder,filesep,'ja_pose_coeff_',ratname,'.mat');
eig_jadyn_coeff_file = strcat(eigenposture_save_folder,filesep,'eig_ja_coeff_',ratname,'.mat');


appearance_coeff_file = strcat(eigenposture_save_folder,filesep,'appearance_coeff_',ratname,'.mat');
dyn_coeff_file_markers = strcat(eigenposture_save_folder,filesep,'dynamics_coeff_2_',ratname,'.mat');
if aligncommonrat
    savefilename = strcat(savefolder,filesep,savename,'_common.mat');
else
    savefilename = strcat(savefolder,filesep,savename,'.mat');
end

if (exist(savefilename,'file') && ~overwrite_savefile)
    fprintf('loading ML features from file %s \n',savefilename)
    tic
    MLmatobj = matfile(savefilename);
    toc
    %load(savefilename)
else
    fprintf('making new ML features from file \n')
    
    
    framelist_true = mocapstruct.modular_cluster_properties.clipped_index{markersetind};%intersect(mocapstruct.modular_cluster_properties.clipped_index{markersetind},framelist);
    
    ML_features.framelist_true = framelist_true;
    %overwrite saved file
    save(savefilename,'-struct','ML_features','-v7.3');
    if numel(framelist_true)>10 %need at least 1 s of data
        
        %% ---------------------------------------------------------------------
        %eigen postures
        agg_pose_features = [];%cat(1,squeeze(marker_velocity(:,1:frame_length_here,4)));
        mocapstruct.modular_cluster_properties.cluster_markersets{markersetind} = [1:10];
        for ll = mocapstruct.modular_cluster_properties.cluster_markersets{markersetind}
            agg_pose_features = cat(1,agg_pose_features,mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{ll})...
                (framelist_true,:)');
        end
        clustering_ind_2 = 1:decimation_factor:size(agg_pose_features,2);
        
        pose_mean = nanmean(agg_pose_features(:, clustering_ind_2),2);
        %center the pose
        centered_pose = bsxfun(@minus,agg_pose_features,pose_mean);
        
        % if doesnt exist compute eigenposture feature coefficients and dynamics
        % coefficients
        if (~exist(pose_coeff_file,'file') || overwrite_coeff)
            [COEFF, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(agg_pose_features');
            save(pose_coeff_file,'COEFF','EXPLAINED');
        else
            load(pose_coeff_file);
        end
        if size(COEFF,1)>30
            COEFF = COEFF(1:30,1:30);
        end
        %compute score
        pose_score = centered_pose'*COEFF;
        pose_score_whitened = pose_score;%bsxfun(@rdivide,pose_score,std(pose_score,[],1));
        ML_features.pose_score = pose_score_whitened(:,1:num_eigenpcs);
        ML_features.pose_coeffs = COEFF(:,1:num_eigenpcs);
        ML_features.pose_mean = pose_mean;
        
        for ll = 1:3
            figure(399)
            h=subplot(1,3,ll)
            visualize_mlfeat_eigenposes(ML_features,mocapstruct,ll,h,markersetind)
        end
        
        
        %% get the local pose features on different timescales --
        timescales = [10,33,100,300];
        for zz = 1:numel(timescales)
            params.gaussorder = timescales(zz)./2;
            gfilter = fspecial('gaussian',[timescales(zz)*6 1], params.gaussorder); %get out to 3 sigma
            ML_features.(strcat('pose_score_',num2str(timescales(zz))))= convn(ML_features.pose_score,gfilter,'same');
        end
        
        
        %open up memory
        clear centered_pose agg_pose_features  pose_score pose_score_whitened
        
        
        %% compute the spectrogram of the whitened eigenfeatures
        %get the clipped time trace of the score
        feature_mat_dyn = [];

        clipped_pre =  hipass_clip_fragments(mocapstruct.markers_aligned_preproc,framelist_true,params);
        maxframes = size(clipped_pre.HeadF,1);
        
        for mm = mocapstruct.modular_cluster_properties.cluster_markersets{markersetind}
            feature_mat_dyn = cat(2,feature_mat_dyn,clipped_pre.(mocapstruct.markernames{mm}));
        end
        feature_mat_dyn_centered = bsxfun(@minus,feature_mat_dyn,mean(feature_mat_dyn,1));
        feature_mat_dyn_centered = feature_mat_dyn_centered*COEFF;        
        
        [dyadic_spectrograms,fr,timespect] = get_dyadic_spectrogram(feature_mat_dyn_centered(:,1:num_eigenpcs)',opts);
        dyadic_spectrograms_reshaped = reshape(dyadic_spectrograms,[],num_eigenpcs,size(dyadic_spectrograms,2));     
        eigenpose_dynamics_coeffs = zeros(num_eigenpcs,numel(fr),num_dynamics_pcs);
        clear feature_mat_dyn_centered
        
        %loop over all spectrograms
        if (~exist(dyn_coeff_file,'file') || overwrite_coeff)
            for ll = 1:size(dyadic_spectrograms_reshaped,2)
                [COEFF_dyn, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(squeeze(dyadic_spectrograms_reshaped(:,ll,:))');
                eigenpose_dynamics_coeffs(ll,:,:) = COEFF_dyn(:,1:num_dynamics_pcs);
            end
            save(dyn_coeff_file,'eigenpose_dynamics_coeffs');
        else
            load(dyn_coeff_file)
        end
        ML_features.eigenpose_dynamics_coeffs = eigenpose_dynamics_coeffs;
        
        
        % compute eigenposture loadings
        dyadic_spectrograms_score = [];
        for ll = 1:size(dyadic_spectrograms_reshaped,2)
            dyadic_spectrograms_white = bsxfun(@minus,squeeze(dyadic_spectrograms_reshaped(:,ll,:)),mean(squeeze(dyadic_spectrograms_reshaped(:,ll,:)),2));
            dyadic_spectrograms_score_temp = dyadic_spectrograms_white'*squeeze( eigenpose_dynamics_coeffs(ll,:,:));
            %dyadic_spectrograms_score_temp = bsxfun(@rdivide,dyadic_spectrograms_score_temp,std(dyadic_spectrograms_score_temp,[],2));
            %plot(dyadic_spectrograms_score_temp)
            dyadic_spectrograms_score = cat(2,dyadic_spectrograms_score,dyadic_spectrograms_score_temp);
        end
        ML_features.dynamics_score = dyadic_spectrograms_score;
        replication_factor = floor(size(feature_mat_dyn,1)./size(dyadic_spectrograms_score,1));
        ML_features.dynamics_replication_ind = reshape(repmat(1:size(dyadic_spectrograms_score,1),replication_factor,1),1,[]);
        
        %open up memory
        clear dyadic_spectrograms_score dyadic_spectrograms_reshaped  dyadic_spectrograms
        
        
        
        
        
        
        %% hand designed pose features
        %high rear
        ML_features.high_rear = mocapstruct.markers_aligned_preproc.SpineF(framelist_true,3)-mocapstruct.markers_aligned_preproc.SpineL(framelist_true,3);
        ML_features.very_high_rear = mocapstruct.markers_aligned_preproc.SpineF(framelist_true,3)-mocapstruct.markers_aligned_preproc.SpineL(framelist_true,3);
        
        %low rear -- shortens stance more and more
        ML_features.low_rear = mocapstruct.markers_aligned_preproc.HeadB(framelist_true,3)-mocapstruct.markers_aligned_preproc.SpineF(framelist_true,3);
        
        %l/r groom
        ML_features.RGroom = mocapstruct.markers_aligned_preproc.SpineF(framelist_true,1)-mocapstruct.markers_aligned_preproc.SpineL(framelist_true,1);
        ML_features.LGroom =mocapstruct.markers_aligned_preproc.SpineL(framelist_true,1)-mocapstruct.markers_aligned_preproc.SpineF(framelist_true,1);
        
        
        % get standard deviation over windows
        
        
        %% morphology features -- inter marker distances
        appearance_pairs = {{'HeadB','SpineF'},{'SpineF','SpineM'},{'SpineM','SpineL'},{'SpineM','Offset1'},{'Offset2','Offset1'},...
            {'Offset2','SpineL'},{'Offset1','SpineF'}};
        
        appearance_features = zeros(size(mocapstruct.markers_preproc.HeadF,1),numel(appearance_pairs));
        markers_to_loop = mocapstruct.modular_cluster_properties.cluster_markersets{markersetind}(1:10);
        appearance_features_agg = zeros(size(mocapstruct.markers_preproc.HeadF(framelist_true,:),1),numel(markers_to_loop )*(numel(markers_to_loop )-1)./2);
        marker_fn =  fieldnames(mocapstruct.markers_preproc);
        
        
        %% Get the PCS of the relative distances
        ind_use = 1;
        for ll = 1:numel(markers_to_loop)
            for jj = ll+1:numel(markers_to_loop)
                
                appearance_features_agg(:,ind_use) = ...
                    vectornorm( mocapstruct.markers_aligned_preproc.(marker_fn{ll})(framelist_true,:),mocapstruct.markers_aligned_preproc.(marker_fn{jj})(framelist_true,:),2);
                ind_use = ind_use +1;
            end
        end
        appearance_features_agg = appearance_features_agg';
        appearance_mean = mean(appearance_features_agg(:,:),2);
        %center the pose
        centered_appearance = bsxfun(@minus,appearance_features_agg,appearance_mean);
        
        if (~exist(appearance_coeff_file,'file') || overwrite_coeff)
            %  for ll = 1:size(dyadic_spectrograms_reshaped,2)
            [COEFF, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(appearance_features_agg(:,:)');
            appearance_coeffs = COEFF;
            % end
            save(appearance_coeff_file,'appearance_coeffs','EXPLAINED');
        else
            load(appearance_coeff_file)
        end
        COEFF = appearance_coeffs;
        
        centered_appearance = centered_appearance'*COEFF;
        ML_features.appearance_features_agg_score_whitened = centered_appearance(:,1:num_eigenpcs);
        
        
        %% get the local pose features on different timescales --
        timescales = [10,33,100,300];
        for zz = 1:numel(timescales)
            params.gaussorder = timescales(zz)./2;
            gfilter = fspecial('gaussian',[timescales(zz)*6 1], params.gaussorder); %get out to 3 sigma
            ML_features.(strcat('appearance_features_agg_score_whitened_',num2str(timescales(zz))))= convn(ML_features.appearance_features_agg_score_whitened,gfilter,'same');
        end
        
        
        clear appearance_features_agg centered_appearance
        
        %% get selected differences
        for ll = 1:numel(appearance_pairs)
            appearance_features(:,ll) = vectornorm( mocapstruct.markers_preproc.(appearance_pairs{ll}{1}),mocapstruct.markers_preproc.(appearance_pairs{ll}{2}),2);
        end
        
        ML_features.appearance_features = appearance_features(framelist_true,:);
        clear appearance_features
        
        
        
        %% save the ML file, clear the features
                if (~overwrite_coeff)
        save(savefilename,'-struct','ML_features','-append','-v7.3');
                end
        ML_features = rmfield(ML_features,fieldnames(ML_features));
        
        
        
        
        
        %% wavelet properties
        downsample = 3;
        frames_use = 1:downsample: maxframes;
        clustering_ind = frames_use; %intersect with the chunking
        cluster_fps = 300./downsample;
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
        opts.params.minF = 0.5; % min freq to analyze
        opts.params.maxF = 60; % max freq to analyze
        
        
        
        %% trajectory features
        % get average trunk displacement, velocity, acceleration, over different
        % time windows: 10,30,100,300 frame windows
        % also head head vel, arm and leg velocity
        clipped_pre_names = fieldnames(clipped_pre);
        
        difforders = [10,33,100,300];
        params.difforder = 10;
        params.medfiltorder = 3;
        params.gaussorder = 2.5;
        ML_features.trunk_vel =zeros(numel(difforders),numel(framelist_true));
        ML_features.head_vel =zeros(numel(difforders),numel(framelist_true));
        
        velcomp_names = {'abs','x','y','z'};
        absolute_velocity_names = {'trunk'};
        absolute_velocity_markers = {[4:8]};
        
        rel_velocity_names = {'head','trunk','hipL','hipR','armL','armR','legL','legR'};
        rel_velocity_markers = {[1:3],[4,6:8],[9],[10],[11,12,13],[14,15,16],[10 17 20],[9 18 19]};
        num_spectrogram_pcs= 15;
        for ll = 1:numel(difforders)
            fprintf('starting absolute and relative velocity for windowsize %f \n',difforders(ll));
            
            params.difforder_movav = difforders(ll);
            for kk = 1:numel(absolute_velocity_names)
                [temp1,velcomp1,velstd,velstd_comp,accel,accel_comp,accelstd,accelstd_comp] = get_markersubgroup_velocity(mocapstruct.markers_preproc,absolute_velocity_markers{kk},params);
                ML_features.(strcat('absolute_velocity_',absolute_velocity_names{kk},'_',velcomp_names{1},'_',num2str(difforders(ll)))) = temp1(framelist_true);
                ML_features.(strcat('absolute_std_velocity_',absolute_velocity_names{kk},'_',velcomp_names{1},'_',num2str(difforders(ll)))) = velstd(framelist_true);
                ML_features.(strcat('absolute_accel_',absolute_velocity_names{kk},'_',velcomp_names{1},'_',num2str(difforders(ll)))) = accel(framelist_true);
                ML_features.(strcat('absolute_std_accel_',absolute_velocity_names{kk},'_',velcomp_names{1},'_',num2str(difforders(ll)))) = accelstd(framelist_true);
                
                for jj =1:size(velcomp1,2)
                    ML_features.(strcat('absolute_velocity_',absolute_velocity_names{kk},'_',velcomp_names{1+jj},'_',num2str(difforders(ll)))) =velcomp1(framelist_true,jj);
                    ML_features.(strcat('absolute_std_velocity',absolute_velocity_names{kk},'_',velcomp_names{1+jj},'_',num2str(difforders(ll)))) =velstd_comp(framelist_true,jj);
                    ML_features.(strcat('absolute_accel',absolute_velocity_names{kk},'_',velcomp_names{1+jj},'_',num2str(difforders(ll)))) =accel_comp(framelist_true,jj);
                    ML_features.(strcat('absolute_std_accel',absolute_velocity_names{kk},'_',velcomp_names{1+jj},'_',num2str(difforders(ll)))) =accelstd_comp(framelist_true,jj);
                    
                end
            end
            
            
            %% relative velocity of the subgroups
            for kk = 1:numel(rel_velocity_names)
                [temp1,velcomp1,velstd,velstd_comp,accel,accel_comp,accelstd,accelstd_comp] = get_markersubgroup_velocity( clipped_pre,rel_velocity_markers{kk},params);
                ML_features.(strcat('rel_velocity_',rel_velocity_names{kk},'_',velcomp_names{1},'_',num2str(difforders(ll)))) = temp1(:);
                ML_features.(strcat('rel_std_velocity_',rel_velocity_names{kk},'_',velcomp_names{1},'_',num2str(difforders(ll)))) = velstd(:);
                
                ML_features.(strcat('rel_accel_',rel_velocity_names{kk},'_',velcomp_names{1},'_',num2str(difforders(ll)))) = accel(:);
                ML_features.(strcat('rel_std_accel_',rel_velocity_names{kk},'_',velcomp_names{1},'_',num2str(difforders(ll)))) = accelstd(:);
                
                for jj =1:size(velcomp1,2)
                    ML_features.(strcat('rel_velocity_',rel_velocity_names{kk},'_',velcomp_names{1+jj},'_',num2str(difforders(ll)))) =velcomp1(:,jj);
                    ML_features.(strcat('rel_std_velocity_',rel_velocity_names{kk},'_',velcomp_names{1+jj},'_',num2str(difforders(ll)))) =velstd_comp(:,jj);
                    ML_features.(strcat('rel_accel_',rel_velocity_names{kk},'_',velcomp_names{1+jj},'_',num2str(difforders(ll)))) =accel_comp(:,jj);
                    ML_features.(strcat('rel_std_accel_',rel_velocity_names{kk},'_',velcomp_names{1+jj},'_',num2str(difforders(ll)))) =accelstd_comp(:,jj);
                end
            end
        end
        
        clear velstd_comp velcomp1
        
        %% spectrograms of the subgroups
        
        %% load old coefficients or not
        if (~exist(dyn_coeff_file_markers,'file') || overwrite_coeff)
            COEFFS_feat = cell(1,numel(rel_velocity_names));
        else
            load(dyn_coeff_file_markers)
        end
        
        %% get the pcs of the spectrogram of markers
        opts.clustering_overlap = 75;
        for kk = 1:4%numel(rel_velocity_names)
            fprintf('computing spectrograms for markers %f \n',kk)
            
            % should refactor this get_dyadic_pcs(clipped_pre,rel_velocity_markers{kk},opts)
            agg_features_here = [];
            for ll = ( rel_velocity_markers{kk})
                agg_features_here = cat(2,agg_features_here,clipped_pre.(clipped_pre_names{ll}));
            end
            %  [dyadic_spectrograms,fr,~] = get_dyadic_spectrogram( agg_features_here',opts);
            dyadic_spectrograms = [];
            params.tapers = [5 7];
            params.Fs = 300;
            

            
            for jj = 1:size(agg_features_here,2)
                if (jj==1)
               [dyadic_spectrograms_temp,fr_temp,tout] = get_dyadic_spectrogram(agg_features_here(:,1)',opts);
               %some weird resizing for small files
               dyadic_spectrograms_temp=dyadic_spectrograms_temp(1:47,:);
               fr_temp = fr_temp(1:47);
               good_freq = find(fr_temp<30);
                else
                
                
                [dyadic_spectrograms_temp,fr_temp,tout] = get_dyadic_spectrogram(agg_features_here(:,jj)',opts);
                               dyadic_spectrograms_temp=dyadic_spectrograms_temp(1:47,:);
                end
                dyadic_spectrograms_temp = dyadic_spectrograms_temp+4;
                dyadic_spectrograms_temp(dyadic_spectrograms_temp<=0) = 0;
                %g_features_here(:,jj),1)-1)./300,agg_features_here(:,jj))
                %     linkaxes(ax,'x')
                
                weighting_function = (5+20*(fr_temp(good_freq ))./30)./5;
                
                dyadic_spectrograms = cat(1,dyadic_spectrograms,bsxfun(@times,dyadic_spectrograms_temp(good_freq,:),weighting_function));
                
                
           end
            % dyadic_spectrograms = log( dyadic_spectrograms);
            
            %% if need new coefficients -- these are constant across files
            if (~exist(dyn_coeff_file_markers,'file') || overwrite_coeff)
                [COEFFS_feat{kk}, dyadic_spectrograms_score, ~, ~,explained] = pca(squeeze(dyadic_spectrograms'));
            else
                if numel(COEFFS_feat{kk})
                dyadic_spectrograms_score = bsxfun(@minus,squeeze(dyadic_spectrograms),mean(squeeze(dyadic_spectrograms),1) )'*squeeze( COEFFS_feat{kk});
                end
            end
            
            if numel(COEFFS_feat{kk})
            replication_factor = floor(size(agg_features_here,1)./size(dyadic_spectrograms_score,1));
            
            %dynamics_pcs = cat(1,dynamics_pcs,zeros(size(agg_features_here,1)-size(dynamics_pcs,1),size(dynamics_pcs,2)));
            
            ML_features.(strcat('spectrogram_pcs_',rel_velocity_names{kk})) = ...
                cat(1,repelem(dyadic_spectrograms_score(:,1:num_spectrogram_pcs),replication_factor,1),...
                zeros(size(agg_features_here,1)-size(dyadic_spectrograms_score,1)*replication_factor,num_spectrogram_pcs));
            ML_features.(strcat('spectrogram_coeffs_',rel_velocity_names{kk})) = COEFFS_feat{kk}(:,1:num_spectrogram_pcs);
            end
        end
        
        if (~exist(dyn_coeff_file_markers,'file') || overwrite_coeff)
            save(dyn_coeff_file_markers,'COEFFS_feat')
        end
        
        figure(109)
        plot3(ML_features.spectrogram_pcs_head(1:300:end,1),ML_features.spectrogram_pcs_head(1:300:end,2),ML_features.spectrogram_pcs_head(1:300:end,3),'+')
        
        
        
        
        %% Get the spectrogram of the velocity
        
        fn_ml = fieldnames(ML_features);
        abs_fn = intersect(intersect(find(cellfun(@numel,strfind(fn_ml,'100'))),find(cellfun(@numel,strfind(fn_ml,'rel_')))),find(cellfun(@numel,strfind(fn_ml,'abs'))));
        abs_fn2 = intersect(find(cellfun(@numel,strfind(fn_ml,'100'))),find(cellfun(@numel,strfind(fn_ml,'_z_'))));
        abs_fn3 = intersect(find(cellfun(@numel,strfind(fn_ml,'100'))),find(cellfun(@numel,strfind(fn_ml,'absolute_velocity_trunk_abs'))));
        
        abs_fn = cat(1,abs_fn,abs_fn2,abs_fn3);
        fn_ml(abs_fn)
        agg_dyn = [];
        for kk = 1:numel(abs_fn)
            agg_dyn = cat(2,agg_dyn,ML_features.(fn_ml{abs_fn(kk)}));
        end
        agg_dyn = agg_dyn';
        dyn_mean = mean(agg_dyn(:,:),2);
        %center the pose
        centered_dyn = bsxfun(@minus,agg_dyn,dyn_mean);
        agg_dyn(abs(agg_dyn)>20) = 20;
        % if doesnt exist compute eigenposture feature coefficients and dynamics
        % % coefficients
        % if (~exist(dyn_coeff_file_2,'file') || overwrite_coeff)
        [COEFF, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(agg_dyn(:,:)');
        %     save(pose_coeff_file,'COEFF');
        % else
        %     load(pose_coeff_file);
        % end
        %compute score
        dyn_score = centered_dyn'*COEFF;
        dyn_score_whitened = dyn_score;%bsxfun(@rdivide,pose_score,std(pose_score,[],1));
        ML_features.dyn_score = dyn_score_whitened(:,1:5);%num_eigenpcs);
        
        figure(44)
        plot(ML_features.dyn_score(:,3))
        %% end euclidean dynamical features
        if (~overwrite_coeff)
       try
           save(savefilename,'-struct','ML_features','-append','-v7.3');
       catch ME
                      save(savefilename,'-struct','ML_features','-append','-v7.3');
       end
        end
        ML_features = rmfield(ML_features,fieldnames(ML_features));
        
        
        
        
        
        
        
        
        
        %% get joint angle features
        %saggital/cross section (ie side view)
        fprintf('computing joing angles \n');
        
        saggital_names = {'head_sagg','neck_sagg','spine_sagg'};
        saggital_pairs =  {[2,3],[3,4],[4,5]}; %head, neck, spine angles , look in the z-y plane
        saggital_include = [1 1 1];
        
        %transverse/overhead
        transverse_names = {'head_trans','neck_trans','spine_trans','hipl_trans','hipr_trans','shouldl_trans','shouldr_trans'};
        transverse_pairs =  {[2,3],[3,4],[4,5],[5,6],[5,7],[4 8 ], [4 9 ],}; %head, neck, spine angles , look in the z-y plane
        transverse_include = [1 1 1 1 1 0 0];
        
        %coronal/along spine (front view)
        coronal_names = {'head_coronal','hipl_coronal','hipr_coronal','shouldl_coronal','shouldr_coronal'};
        coronal_pairs =  {[1,3],[5,6],[5,7],[4 8], [4 9]}; %head, neck, spine angles , look in the z-y plane
        coronal_include = [1 1 1 0 0];
        
        %alljt names
        allangles_names = {'lelbow_all','larm_all','relbow_all','rarm_all','lknee_all','lshin_all','rknee_all','rshin_all'};
        allangles_pairs =  {[8,10],[10,11],[9,12],[12,13],[6,14],[14,15],[7,16],[16,17]}; %head, neck, spine angles , look in the z-y plane
        all_include = zeros(1,numel(allangles_names));
        
        %% specify the specific angles for the different appendages
        appendage_names = {'Head','LArm','Rarm','trunk','LLeg','RLeg'};
        appendage_anglegps = cell(1,numel(appendage_names));
        appendage_anglegps{1} = {'head_sagg','neck_sagg','head_trans','neck_trans','head_coronal'};
        appendage_anglegps{2} = {'shouldl_trans','shouldl_coronal','lelbow_all','larm_all'};
        appendage_anglegps{3} = {'shouldr_trans','shouldr_coronal','relbow_all','rarm_all'};
        appendage_anglegps{4} = {'spine_sagg','spine_trans'};
        appendage_anglegps{5} = {'hipl_trans','hipl_coronal','lknee_all','lshin_all'};
        appendage_anglegps{6} = {'hipr_trans','hipr_coronal','rknee_all','rshin_all'};
        
        ML_features.appendage_names = appendage_names;
        
        %% appendage segment lengths
        appendage_segvals = cell(1,numel(appendage_names));
        appendage_segvals{1} = [1,2,3];
        appendage_segvals{2} = [8,10,11];
        appendage_segvals{3} = [9,12,13];%'shouldr_trans','shouldr_coronal','relbow_all','rarm_all'};
        appendage_segvals{4} = [4,5];
        appendage_segvals{5} = [6,14,15];%'hipl_trans','hipl_coronal','lknee_all','lshin_all'};
        appendage_segvals{6} = [7,16,17];
        
        
        saggital_inds = [2,3];
        coronal_inds = [1,3];
        transverse_inds = [1,2];
        allangles_inds = [1,2,3]; %use on knees and arms
        %transverse_pairs
        
        %% get the various
        segment_pairs = {{'HeadB','HeadL'},{'HeadF','HeadB'},{'HeadB','SpineF'},{'SpineF','SpineM'} ,...%1-4
            {'SpineL','SpineM'},{'SpineL','HipL'},{'SpineL','HipR'},... %5-7
            {'SpineF','ShoulderL'},{'SpineF','ShoulderR'},... %8,9
            {'ShoulderL','ElbowL'},{'ElbowL','ArmL'},{'ShoulderR','ElbowR'},{'ElbowR','ArmR'},...%10-13
            {'HipL','KneeL'},{'KneeL','ShinL'},{'HipR','KneeR'},{'KneeR','ShinR'}}; %14-17
        
        jointangle_struct = struct();
        mean_seglengths = struct();
        all_seglengths = cell(1,numel(segment_pairs));
        
        
        for ll = 1:numel(saggital_pairs)
            vec1 =  mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(1)}{1})-...
                mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(1)}{2});
            
            vec2 =  mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(2)}{1})-...
                mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(2)}{2});
            % jointangle_struct.saggital_names{ll} = [];
            % if saggital_include(ll)
            jointangle_struct.(saggital_names{ll}) = ...
                acosd(dot(vec1(framelist_true,saggital_inds)',vec2(framelist_true,saggital_inds)')...
                ./(sqrt(sum(vec1(framelist_true,saggital_inds).^2,2)).*sqrt(sum(vec2(framelist_true,saggital_inds).^2,2)))')';
            mean_seglengths.(saggital_names{ll}) = [nanmean(sqrt(sum(vec1(framelist_true,saggital_inds).^2,2))) ...
                nanmean(sqrt(sum(vec2(framelist_true,saggital_inds).^2,2)))];
            % all_seglengths.(saggital_names{ll}) = [(sqrt(sum(vec1(framelist_true,saggital_inds).^2,2))) ...
            %   (sqrt(sum(vec2(framelist_true,saggital_inds).^2,2)))];
            %  end
            
        end
        
        for ll = 1:numel(coronal_pairs)
            vec1 =  mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(1)}{1})-...
                mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(1)}{2});
            
            vec2 =  mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(2)}{1})-...
                mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(2)}{2});
            % jointangle_struct.saggital_names{ll} = [];
            %  if coronal_include(ll)
            
            jointangle_struct.(coronal_names{ll}) = ...
                acosd(dot(vec1(framelist_true,coronal_inds)',vec2(framelist_true,coronal_inds)')...
                ./(sqrt(sum(vec1(framelist_true,coronal_inds).^2,2)).*sqrt(sum(vec2(framelist_true,coronal_inds).^2,2)))')';
            mean_seglengths.(coronal_names{ll}) = [nanmean(sqrt(sum(vec1(framelist_true,coronal_inds).^2,2))) ...
                nanmean(sqrt(sum(vec2(framelist_true,coronal_inds).^2,2)))];
            %  end
            
        end
        
        
          
        %% all angles
        for ll = 1:numel(allangles_pairs)
            vec1 =  mocapstruct.markers_aligned_preproc.(segment_pairs{allangles_pairs{ll}(1)}{1})-...
                mocapstruct.markers_aligned_preproc.(segment_pairs{allangles_pairs{ll}(1)}{2});
            
            vec2 =  mocapstruct.markers_aligned_preproc.(segment_pairs{allangles_pairs{ll}(2)}{1})-...
                mocapstruct.markers_aligned_preproc.(segment_pairs{allangles_pairs{ll}(2)}{2});
            % jointangle_struct.saggital_names{ll} = [];
            %  if all_include(ll)
            jointangle_struct.(allangles_names{ll}) = ...
                acosd(dot(vec1(framelist_true,allangles_inds)',vec2(framelist_true,allangles_inds)')...
                ./(sqrt(sum(vec1(framelist_true,allangles_inds).^2,2)).*sqrt(sum(vec2(framelist_true,allangles_inds).^2,2)))')';
            mean_seglengths.(allangles_names{ll}) = [nanmean(sqrt(sum(vec1(framelist_true,allangles_inds).^2,2))) ...
                nanmean(sqrt(sum(vec2(framelist_true,allangles_inds).^2,2)))];
            % end
            
        end
        
        %% transverse
        for ll = 1:numel(transverse_pairs)
            vec1 =  mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(1)}{1})-...
                mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(1)}{2});
            
            vec2 =  mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(2)}{1})-...
                mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(2)}{2});
            % jointangle_struct.saggital_names{ll} = [];
            %if transverse_include(ll)
            jointangle_struct.(transverse_names{ll}) = ...
                acosd(dot(vec1(framelist_true,transverse_inds)',vec2(framelist_true,transverse_inds)')...
                ./(sqrt(sum(vec1(framelist_true,transverse_inds).^2,2)).*sqrt(sum(vec2(framelist_true,transverse_inds).^2,2)))')';
            mean_seglengths.(transverse_names{ll}) = real([nanmean(sqrt(sum(vec1(framelist_true,transverse_inds).^2,2))) ...
                nanmean(sqrt(sum(vec2(framelist_true,transverse_inds).^2,2)))]);
            % end
            
        end
      
        
        %% get the segment lengths
        for ll = 1:numel(segment_pairs)
            vec1 =  mocapstruct.markers_aligned_preproc.(segment_pairs{ll}{1})-...
                mocapstruct.markers_aligned_preproc.(segment_pairs{ll}{2});
            
            all_seglengths{ll} = (sqrt(sum(vec1(framelist_true,allangles_inds).^2,2)));
            % end
        end
        
        
        total_include = cat(2,saggital_include,coronal_include,all_include,transverse_include);
        fn_ja_all = fieldnames(jointangle_struct);
        
        for ll =1:numel(fn_ja_all)
             jointangle_struct.(fn_ja_all{ll})(isnan(jointangle_struct.(fn_ja_all{ll}))) = 0;
             jointangle_struct.(fn_ja_all{ll})(isinf(jointangle_struct.(fn_ja_all{ll}))) = 0;
        end

        
        jointangle_struct = structfun(@(x) real(x),jointangle_struct, 'UniformOutput', false);
        mean_seglengths = structfun(@(x) real(x),mean_seglengths, 'UniformOutput', false);
        
        ML_features.joint_angles_mean = real(structfun(@nanmean,(jointangle_struct)));
        
        
        
        %% hipass and clip all of the joint angles
        jointangle_struct_hipass=hipass_clip_fragments(jointangle_struct,1:numel(jointangle_struct.head_sagg),params);
        
        params.difforder = 10;
        params.medfiltorder = 3;
        params.gaussorder = 2.5;
        jointangle_struct_hipass= smooth_struct(jointangle_struct_hipass,params);
        ML_features.joint_angles = jointangle_struct_hipass;
        
        ML_features.mean_seg_lengths = reshape(struct2array(mean_seglengths),[numel(fieldnames(ML_features.joint_angles)) 2]);
        
        figure(66)
        plot(( real(jointangle_struct.larm_all)))
        
        %jointangle_struct_hipass.(fn_ja_all)
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
        
        for kk = 1:numel(appendage_anglegps)
            appendage_anglevals{kk} = [];
            %%get the angles and run PCA
            for ll = 1:numel(appendage_anglegps{kk})
                index = find(strcmp( fn_ja_all ,appendage_anglegps{kk}{ll}));
                appendage_anglevals{kk} = cat(2,appendage_anglevals{kk},jointangle_struct.(fn_ja_all{index}));
            end
            
            %subtract mean
            appendage_anglevals{kk} = bsxfun(@minus,appendage_anglevals{kk},nanmean(appendage_anglevals{kk},1));
            [COEFFS_appendages{kk}, ~, ~, ~,appendage_explained{kk}] = pca(squeeze(appendage_anglevals{kk}));
            %   save(ja_pose_coeffs_file,'COEFFS','explained');
            if numel(COEFFS_appendages{kk})
                appendage_dyadic_spectrograms{kk} =  appendage_anglevals{kk}*COEFFS_appendages{kk};
            end
            
            %% do for segments
            appendage_lengths{kk} = [];
            for ll =1:numel(appendage_segvals{kk})
                appendage_lengths{kk} = cat(2,appendage_lengths{kk},all_seglengths{appendage_segvals{kk}(ll)});
            end
            
              %subtract mean
            appendage_lengths{kk} = bsxfun(@minus,appendage_lengths{kk},nanmean(appendage_lengths{kk},1));
            [COEFFS_appendages_lengths{kk}, ~, ~, ~,appendage_explained_lengths{kk}] = pca(squeeze(appendage_lengths{kk}));
            %   save(ja_pose_coeffs_file,'COEFFS','explained');
            if numel(COEFFS_appendages_lengths{kk})
                appendage_dyadic_spectrograms_lengths{kk} =  appendage_lengths{kk}*COEFFS_appendages_lengths{kk};
            end
            
        end
        
        
        ML_features.appendage_coeffs = COEFFS_appendages;
        ML_features.appendage_pca_score = appendage_dyadic_spectrograms;
        ML_features.appendage_pca_explained = appendage_explained;
        
          
        
        ML_features.appendage_coeffs_lengths = COEFFS_appendages_lengths;
        ML_features.appendage_pca_score_lengths = appendage_dyadic_spectrograms_lengths;
        ML_features.appendage_pca_explained_lengths = appendage_explained_lengths;
        
        
        %% get the appendage segment length pcs as well
        
        %% pcs of joint angles of total top set
        principal_jas = total_include;
        fn_ja = fn_ja_all(find(principal_jas));
        agg_ja_features = [];
        for ll = 1:numel(fn_ja)
            agg_ja_features = cat(2,agg_ja_features,jointangle_struct.(fn_ja{ll}));
        end
        agg_ja_features = bsxfun(@minus,agg_ja_features,nanmean(agg_ja_features,1));
        agg_ja_features(isnan(agg_ja_features)) = 0;
        agg_ja_features(isinf(agg_ja_features)) = 0;
        
        if (~exist(ja_pose_coeffs_file,'file') || overwrite_coeff)
            [COEFFS, dyadic_spectrograms, ~, ~,explained] = pca(squeeze(agg_ja_features));
            save(ja_pose_coeffs_file,'COEFFS','explained');
        else
            load(ja_pose_coeffs_file);
        end
        dyadic_spectrograms = agg_ja_features*COEFFS;
        
        ML_features.ja_pca_coeffs = COEFF;
        ML_features.ja_dyadic_spectrograms = dyadic_spectrograms(:,1:end);
        ML_features.ja_pca_explained = explained;
        
        %% smooth the dyadic spectrograms as well -- get the pose over a local region
        timescales = [10,33,100,300];
        for zz = 1:numel(timescales)
            params.gaussorder = timescales(zz)./2;
            gfilter = fspecial('gaussian',[timescales(zz)*6 1], params.gaussorder); %get out to 3 sigma
            ML_features.(strcat('ja_dyadic_spectrograms_',num2str(timescales(zz))))= convn(ML_features.ja_dyadic_spectrograms,gfilter,'same');
            
            appendage_pca_score_smoothed = cell(1,numel(appendage_anglegps));
            for kk = 1:numel(appendage_anglegps)
                appendage_pca_score_smoothed{kk} = convn(ML_features.appendage_pca_score{kk},gfilter,'same');
            end
            ML_features.(strcat('appendage_pca_score',num2str(timescales(zz)))) = appendage_pca_score_smoothed;
        end
        
        
        
        
        %% load old coefficients or not
        if (~exist(eig_jadyn_coeff_file,'file') || overwrite_coeff)
            COEFFS_feat = cell(1,1);
            COEFFS_feat_wl= cell(1,1);
            explained = cell(1,1);
            explained_wl= cell(1,1);
        else
            load(eig_jadyn_coeff_file)
        end
        
        
        %% look at pcs of the eigenangles
        agg_features=[];
        agg_features_wl=[];
        
        %% JA dyadic spectrogram wavelet and spectorgram coefficients
        for ll = 1:size(dyadic_spectrograms,2)
            dHipass = designfilt('highpassiir', 'FilterOrder', 3, 'HalfPowerFrequency', 0.5/(params.fps/2), ...
                'DesignMethod', 'butter');
            [f1_hipass,f2_hipass] = tf(dHipass);
            frame_fragment = filtfilt(f1_hipass,f2_hipass,dyadic_spectrograms(:,ll));
            
            %% compute spectrograms
            %figure(101)
            % plot(frame_fragment)
            %     [S,t,f]=mtspecgramc(frame_fragment,[1 0.5],params);
            %     dyadic_spectrograms_temp = (S(:,1:50).^2)';
            
            [dyadic_spectrograms_temp,fr_temp,tout] = get_dyadic_spectrogram(frame_fragment(:,1)',opts);
            zsc_dy = exp(dyadic_spectrograms_temp(:,:));
            zsc_dy = bsxfun(@rdivide,bsxfun(@minus,zsc_dy,median(zsc_dy,2)),mad(zsc_dy,[],2));%,mad(zsc_dy,[],2));%zsc_dy
            % zsc_dy = exp(dyadic_spectrograms_temp(:,:));
            % zsc_dy(zsc_dy>1) = 1;
            %zsc_dy = bsxfun(@rdivide,zsc_dy,sum(zsc_dy,1));
            zsc_dy(zsc_dy>5) = 5;
            zsc_dy(zsc_dy<-1) = -1;
            good_freq = find(fr_temp<30);
            
            agg_features = cat(1,agg_features,zsc_dy(good_freq,:));
            
            frame_fragment(find(isnan(frame_fragment))) = 0;
            frame_fragment(find(isinf(frame_fragment))) = 0;
            %% compute the wavelets coefficients
            [wavelet_coeffs,w_map,fr_wavelet]= return_wavelets(frame_fragment(1:3:end,1),1:size(frame_fragment(1:3:end,1),1),opts);
            
            %need to motivate the thresholding
            w_map = w_map+3;
            w_map(w_map<0) = 0;
            agg_features_wl = cat(2,agg_features_wl,w_map);
            
            %% also compute for each of the appendages
        end
        
        
        
        %% compute the PCS of the wavelet
        if (~exist(eig_jadyn_coeff_file,'file') || overwrite_coeff)
            [COEFFS_feat{1}, dyadic_spectrograms_score, ~, ~,explained] = pca(agg_features');
            [COEFFS_feat_wl{1}, dyadic_spectrograms_score_wl, ~, ~,explained_wl{1}] = pca(squeeze(agg_features_wl));
            %dyadic_spectrograms_score = bsxfun(@minus,squeeze(dyadic_spectrograms)',mean(squeeze(dyadic_spectrograms)',2) )*squeeze( COEFFS_feat{kk});
        else
            dyadic_spectrograms_score = bsxfun(@minus,squeeze(agg_features),mean(squeeze(agg_features),2) )'*squeeze( COEFFS_feat{1});
            dyadic_spectrograms_score_wl = bsxfun(@minus,squeeze(agg_features_wl),mean(squeeze(agg_features_wl),2) )*squeeze( COEFFS_feat_wl{1});
        end
        
        %% expand these back out
        replication_factor = floor(size(frame_fragment,1)./size(dyadic_spectrograms_score,1));
        dynamics_pcs = repelem(dyadic_spectrograms_score(:,1:num_spectrogram_pcs),replication_factor,1);
        dynamics_pcs = cat(1,dynamics_pcs,zeros(size(frame_fragment,1)-size(dynamics_pcs,1),size(dynamics_pcs,2)));
        
        %replicate the wavelet
        replication_factor_wl = ceil(size(frame_fragment,1)./size(dyadic_spectrograms_score_wl,1));
        dynamics_pcs_wl = repelem(dyadic_spectrograms_score_wl(:,1:num_spectrogram_pcs),replication_factor_wl,1);
        if (size(dynamics_pcs_wl,1)<=size(frame_fragment,1))
            dynamics_pcs_wl = cat(1,dynamics_pcs_wl,zeros(size(frame_fragment,1)-size(dynamics_pcs_wl,1),size(dynamics_pcs_wl,2)));
        else
            dynamics_pcs_wl((end-(size(dynamics_pcs_wl,1)-size(frame_fragment,1) )):end,:) = [];
        end
        
        
        ML_features.ja_eig_spectrogrampcs = dynamics_pcs(:,1:15);
        ML_features.ja_eig_wlpcs = dynamics_pcs_wl(:,1:15);
        
        ML_features.ja_eig_spectrogrampcs_explained = explained;
        ML_features.ja_eig_wlpcs_explained = explained_wl;
        
        ML_features.ja_eig_spectrogrampcs_coeff = COEFFS_feat;
        ML_features.ja_eig_wlpcs_coeff = COEFFS_feat_wl;
        
        %       zvals =  tsne(  cat(2, ML_features.ja_dyadic_spectrograms(1:25:end,:), ML_features.ja_eig_wlpcs(1:25:end,:)));
        %         figure(66)
        %plot(zvals(:,1),zvals(:,2),'+')
        
        
        %% also compute the wavelet for the appendages
        COEFFS_feat_wl_appendages = cell(1,numel(appendage_anglegps));
        dyadic_spectrograms_score_wl_appendages= cell(1,numel(appendage_anglegps));
        explained_wl_appendages = cell(1,numel(appendage_anglegps));
        
        for kk = 1:numel(appendage_anglegps)
            agg_features_wl = [];
            for ll = 1:size(ML_features.appendage_pca_score{kk},2)
                ML_features.appendage_pca_score{kk}(find(isnan(ML_features.appendage_pca_score{kk}(:,ll))),ll) = 0;
                ML_features.appendage_pca_score{kk}(find(isinf(ML_features.appendage_pca_score{kk}(:,ll))),ll) = 0;
                
                frame_fragment = filtfilt(f1_hipass,f2_hipass,ML_features.appendage_pca_score{kk}(:,ll));
                frame_fragment(find(isnan(frame_fragment))) =0;
                frame_fragment(find(isinf(frame_fragment))) =0;
                
                [wavelet_coeffs,w_map,fr_wavelet]= return_wavelets(frame_fragment(1:3:end,1),1:size(frame_fragment(1:3:end,1),1),opts);
                w_map = w_map+3;
                w_map(w_map<0) = 0;
                agg_features_wl = cat(2,agg_features_wl,w_map);
            end
            [COEFFS_feat_wl_appendages{kk}, dyadic_spectrograms_score_wl_appendages{kk}, ~, ~,explained_wl_appendages{kk}] = pca(squeeze(agg_features_wl));
            
            % replicate elements to fill in
            replication_factor_wl = ceil(size(frame_fragment,1)./size(dyadic_spectrograms_score_wl_appendages{kk}(:,1),1));
            num_spectrogram_pcs =min(num_spectrogram_pcs,size(dyadic_spectrograms_score_wl_appendages{kk},2));
            dyadic_spectrograms_score_wl_appendages{kk} = repelem( dyadic_spectrograms_score_wl_appendages{kk}(:,1:num_spectrogram_pcs),replication_factor_wl,1);
            
            %do a last pruning
            if (size(dyadic_spectrograms_score_wl_appendages{kk},1)<size(frame_fragment,1))
                dyadic_spectrograms_score_wl_appendages{kk} = cat(1,dyadic_spectrograms_score_wl_appendages{kk},zeros(size(frame_fragment,1)...
                    -size(dyadic_spectrograms_score_wl_appendages{kk},1),size(dyadic_spectrograms_score_wl_appendages{kk},2)));
            else
                dyadic_spectrograms_score_wl_appendages{kk}((end-(size(dyadic_spectrograms_score_wl_appendages{kk},1)-size(frame_fragment,1) )):end,:) = [];
            end
        end
        
        ML_features.COEFFS_feat_wl_appendages = COEFFS_feat_wl_appendages;
        ML_features.dyadic_spectrograms_score_wl_appendages = dyadic_spectrograms_score_wl_appendages;
        ML_features.explained_wl_appendages = explained_wl_appendages;
        
        %% save coefficients
        if (~exist( eig_jadyn_coeff_file,'file') || overwrite_coeff)
            save( eig_jadyn_coeff_file,'COEFFS_feat','explained','COEFFS_feat_wl','explained_wl')
            
        end
        
        
        
        
        %% dynamics of the joint angles -- are they cleaner?
        head_angles = {'head_sagg','head_trans','neck_sagg','neck_trans','head_coronal'};
        trunk_angles = {'spine_sagg','hipl_coronal','hipr_coronal','spine_trans','hipl_trans','hipr_trans'};
        
        angle_lists = {head_angles,trunk_angles};
        angle_list_name = {'head_angle','trunk_angle'};
        ML_features.trunk_angles = trunk_angles;
        
        
        
        
        %% load old coefficients or not for indiv angles
        if (~exist(jadyn_coeff_file,'file') || overwrite_coeff)
            COEFFS_feat = cell(1,numel(angle_lists));
            COEFFS_feat_wl= cell(1,numel(angle_lists));
            explained = cell(1,numel(angle_lists));
            explained_wl= cell(1,numel(angle_lists));
        else
            load(jadyn_coeff_file)
        end
        
        
        %% Compute the JA velocity
        difforders = [100,300];
        params.difforder = 10;
        params.medfiltorder = 3;
        params.gaussorder = 2.5;
        for kk = 1:numel(angle_lists)
            for zz = 1:numel(difforders)
                
                ja_velocity_here = [];
                ja_std_here = [];
                ja_accel = [];
                ja_accel_std = [];
                
                for ll = 1:numel( angle_lists{kk})
                    params.difforder_movav = difforders(zz);
                    
                    [temp1,~,velstd,~,accel,~,accelstd,~] = get_markersubgroup_velocity( ML_features.joint_angles,find(strcmp(fieldnames(ML_features.joint_angles) ,...
                        (angle_lists{kk}{ll}))),params);
                    ja_velocity_here = cat(2,  ja_velocity_here ,temp1);
                    ja_std_here = cat(2,  ja_std_here ,velstd);
                    ja_accel = cat(2,  ja_accel ,accel);
                    ja_accel_std = cat(2,  ja_accel_std ,accelstd);
                    
                end
                ML_features.(strcat('ja_velocity',angle_list_name{kk},'_difforder_',num2str(difforders(zz)))) = ja_velocity_here;
                ML_features.(strcat('ja_std',angle_list_name{kk},'_difforder_',num2str(difforders(zz)))) = ja_std_here;
                ML_features.(strcat('ja_accel',angle_list_name{kk},'_difforder_',num2str(difforders(zz)))) = ja_accel;
                ML_features.(strcat('ja_accel_std',angle_list_name{kk},'_difforder_',num2str(difforders(zz)))) = ja_accel_std;
            end
        end
        
        %% compute the marker spectrograms
        opts.clustering_overlap = 75;
        for kk = 1:numel(angle_lists)
            ML_features.angle_names{kk} = (angle_lists{kk});
            fprintf('computing spectrograms for markers %f \n',kk)
            
            % should refactor this get_dyadic_pcs(clipped_pre,rel_velocity_markers{kk},opts)
            agg_features_here = [];
            dyadic_spectrograms = [];
            dyadic_spectrograms_wl = [];
            
            for ll = 1:numel( angle_lists{kk})
                
                
                spectrogram_thresh = 5;
                [dyadic_spectrograms_temp,fr_temp,tout] = get_dyadic_spectrogram(ML_features.joint_angles.(angle_lists{kk}{ll})',opts);
                dyadic_spectrograms_temp = dyadic_spectrograms_temp+spectrogram_thresh;
                dyadic_spectrograms_temp(dyadic_spectrograms_temp<=0) = 0;
                
                good_freq = find(fr_temp<30);
                weighting_function = (5+20*(fr_temp(good_freq ))./30)./5;
                zsc_dy = dyadic_spectrograms_temp;
                %   zsc_dy = bsxfun(@rdivide,zsc_dy,sum(zsc_dy,1));
                
                
                dyadic_spectrograms = cat(1,dyadic_spectrograms,bsxfun(@times,dyadic_spectrograms_temp(good_freq,:),weighting_function));
                
                
                %% compute the wavelets coefficients
                [wavelet_coeffs,w_map,fr_wavelet]= return_wavelets(cat(2,ML_features.joint_angles.(angle_lists{kk}{ll})(1:3:end)),...
                    1:numel(1:3:numel(ML_features.joint_angles.(angle_lists{kk}{ll}))),opts);
                
                w_map = w_map+3;
                w_map(w_map<0) = 0;
                
                %[COEFFS_feat{kk}, dyadic_spectrograms_score, ~, ~,explained] = pca(squeeze(w_map));
                dyadic_spectrograms_wl = cat(2,dyadic_spectrograms_wl,w_map);
                
                
                
            end
            
            %% if need new coefficients -- these are constant across files
            if (~exist(jadyn_coeff_file,'file') || overwrite_coeff)
                [COEFFS_feat{kk}, dyadic_spectrograms_score, ~, ~,explained{kk}] = pca(squeeze(dyadic_spectrograms)');
                [COEFFS_feat_wl{kk}, dyadic_spectrograms_score_wl, ~, ~,explained_wl{kk}] = pca(squeeze(dyadic_spectrograms_wl));
                %dyadic_spectrograms_score = bsxfun(@minus,squeeze(dyadic_spectrograms)',mean(squeeze(dyadic_spectrograms)',2) )*squeeze( COEFFS_feat{kk});
            else
                
                dyadic_spectrograms_score = bsxfun(@minus,squeeze(dyadic_spectrograms),mean(squeeze(dyadic_spectrograms),2) )'*squeeze( COEFFS_feat{kk});
                dyadic_spectrograms_score_wl = bsxfun(@minus,squeeze(dyadic_spectrograms_wl),mean(squeeze(dyadic_spectrograms_wl),2) )*squeeze( COEFFS_feat_wl{kk});
                
            end
            
            replication_factor = floor(size(ML_features.joint_angles.(angle_lists{kk}{ll}),1)./size(dyadic_spectrograms_score,1));
            dynamics_pcs = repelem(dyadic_spectrograms_score(:,1:num_spectrogram_pcs),replication_factor,1);
            dynamics_pcs = cat(1,dynamics_pcs,zeros(size(ML_features.joint_angles.(angle_lists{kk}{ll}),1)-size(dynamics_pcs,1),size(dynamics_pcs,2)));
            
            
            %replicate the wavelet
            replication_factor_wl = ceil(size(ML_features.joint_angles.(angle_lists{kk}{ll}),1)./size(dyadic_spectrograms_score_wl,1));
            dynamics_pcs_wl = repelem(dyadic_spectrograms_score_wl(:,1:num_spectrogram_pcs),replication_factor_wl,1);
            if (size(dynamics_pcs_wl,1)<size(ML_features.joint_angles.(angle_lists{kk}{ll}),1))
                dynamics_pcs_wl = cat(1,dynamics_pcs_wl,zeros(size(ML_features.joint_angles.(angle_lists{kk}{ll}),1)-size(dynamics_pcs_wl,1),size(dynamics_pcs_wl,2)));
            else
                dynamics_pcs_wl((end-(size(dynamics_pcs_wl,1)-size(ML_features.joint_angles.(angle_lists{kk}{ll}),1) )):end,:) = [];
            end
            
            num_spectrogram_pcs = 25;
            ML_features.ja_freq = fr_temp(good_freq);
            ML_features.(strcat('spectrogram_pcs_',angle_list_name{kk},'_explained')) = explained{kk};
            ML_features.(strcat('spectrogram_pcs_wl_',angle_list_name{kk},'_explained_wl')) = explained_wl{kk};
            
            ML_features.(strcat('spectrogram_pcs_',angle_list_name{kk})) = dynamics_pcs;
            ML_features.(strcat('spectrogram_pcs_wl_',angle_list_name{kk})) = dynamics_pcs_wl;
            ML_features.fr_wavelet = fr_wavelet;
            ML_features.(strcat('wavelet_coeffs_',angle_list_name{kk})) = COEFFS_feat_wl{kk}(:,1:num_spectrogram_pcs);
            
            
            ML_features.(strcat('spectrogram_coeffs_',angle_list_name{kk})) = COEFFS_feat{kk}(:,1:num_spectrogram_pcs);
            ML_features.ja_dyn_explained = explained{kk};
        end
        
        if (~exist( jadyn_coeff_file,'file') || overwrite_coeff)
            save( jadyn_coeff_file,'COEFFS_feat','explained','COEFFS_feat_wl','explained_wl')
        end
        
        
        
        %% visualize the spectrogram and wavelet cofficients
        num_pc = 6;
        num_angle = 5;
        
        for pc_plot = 1:6
            kk_plot =2;
            
            COEFFS_resh= reshape(COEFFS_feat{kk_plot}(:,pc_plot),numel( ML_features.ja_freq),[]);
            deviation = nanstd(ML_features.(strcat('spectrogram_pcs_',angle_list_name{kk_plot}))(:,pc_plot),[],1);
            
            summed_coeffs =  COEFFS_resh*deviation;%ML_features.mean_ja_spect{kk_plot}'+COEFFS_resh*deviation;
            summed_coeffs_minus =  -COEFFS_resh*deviation;%ML_features.mean_ja_spect{kk_plot}'-COEFFS_resh*deviation;
            
            summed_coeffs_exp = summed_coeffs;% (10.^(ML_features.mean_ja_spect{kk_plot}'+COEFFS_resh*deviation)-10.^(ML_features.mean_ja_spect{kk_plot}'));
            
            for angle_plot = 1:5
                figure(44)
                subplot(num_angle,num_pc,pc_plot+num_pc*(angle_plot-1))
                plot(ML_features.ja_freq,COEFFS_resh(:,angle_plot)*deviation,'r');
                hold on
                %plot(ML_features.mean_ja_spect{kk_plot}(angle_plot,:)','k');
                plot(ML_features.ja_freq,-COEFFS_resh(:,angle_plot)*deviation,'b' )
                
                if (pc_plot == 1)
                    ylabel(ML_features.angle_names{kk_plot}{angle_plot})
                end
                if (angle_plot == 1)
                    ntitle(strcat('PC ',num2str(pc_plot)));
                end
                
                
                
                figure(45)
                subplot(num_angle,num_pc,pc_plot+num_pc*(angle_plot-1))
                plot(ML_features.ja_freq,summed_coeffs(:,angle_plot),'r');
                hold on
                % plot(ML_features.mean_ja_spect{kk_plot}(angle_plot,:)','k');
                plot(ML_features.ja_freq,summed_coeffs_minus(:,angle_plot),'b' )
                
                if (pc_plot == 1)
                    ylabel(ML_features.angle_names{kk_plot}{angle_plot})
                end
                if (angle_plot == 1)
                    ntitle(strcat('PC ',num2str(pc_plot)));
                end
                
                figure(46)
                subplot(num_angle,num_pc,pc_plot+num_pc*(angle_plot-1))
                plot(ML_features.ja_freq,exp(summed_coeffs(:,angle_plot)),'r');
                hold on
                % plot(exp(ML_features.mean_ja_spect{kk_plot}(angle_plot,:)'),'k');
                plot(ML_features.ja_freq,exp(summed_coeffs_minus(:,angle_plot)),'b');
                
                if (pc_plot == 1)
                    ylabel(ML_features.angle_names{kk_plot}{angle_plot})
                end
                if (angle_plot == 1)
                    ntitle(strcat('PC ',num2str(pc_plot)));
                end
            end
        end
        
        
        
        %
        %
        figure(98)
        %plot3(ML_features.spectrogram_pcs_trunk_angle(1:100:100*10000,1),ML_features.spectrogram_pcs_trunk_angle(1:100:100*10000,2),ML_features.spectrogram_pcs_trunk_angle(1:100:100*10000,3),'+r')
        
        % figure(108)
        %  figure(108)
        %  plot3(ML_features.spectrogram_pcs_trunk_angle(1:100:100*10000,1),ML_features.spectrogram_pcs_trunk_angle(1:100:100*10000,2),ML_features.spectrogram_pcs_trunk_angle(1:100:100*10000,3),'+r')
        %
        figure(99)
        
        plot3(dyadic_spectrograms_score_wl(1:100:end,1),dyadic_spectrograms_score_wl(1:100:end,2),dyadic_spectrograms_score_wl(1:100:end,3),'+')
        
        MLmatobj = matfile(savefilename);
        
        %% save after JA features
                if (~overwrite_coeff)
        save(savefilename,'-struct','ML_features','-append','-v7.3');
                end
        ML_features = rmfield(ML_features,fieldnames(ML_features));
        
        
        
        
        %% create mocap features for knees and arms -- get the eigenposes (relative to a base marker), eigen JA, dynamics and velocity characteristics (for the table)
        %mocapstruct.modular_cluster_properties.cluster_markersets{3}
        % get indicies for each set
        
        % get eigenpos
        
        
        
        
    else
        fprintf('no frames present!!! \n')
        MLmatobj = matfile(savefilename);
        
    end
end
end



%% code for visualization

%xrec = icwt(w_map,300,fr_wavelet,'wname','Mortlet') ;

%    figure(100)
%     ax(1)= subplot(2,1,1)
%   %  imagesc(tout,fr_temp(good_freq),bsxfun(@times,dyadic_spectrograms_temp(good_freq,:),weighting_function));%bsxfun(@rdivide,zsc_dy,sum(zsc_dy,1)))
%       % imagesc(tout,fr_temp(good_freq),bsxfun(@times,zsc_dy(good_freq,:),1));%bsxfun(@rdivide,zsc_dy,sum(zsc_dy,1)))
%     imagesc(0:1./100:(size(ML_features.joint_angles.(angle_lists{kk}{ll})',2)-1)./300,1:size(w_map,2),w_map')
%
%     % caxis([0 0.1])
%     ax(2)= subplot(2,1,2)
%     plot(0:1./300:(size(ML_features.joint_angles.(angle_lists{kk}{ll})',2)-1)./300,ML_features.joint_angles.(angle_lists{kk}{ll}))
%     linkaxes(ax,'x')
% %
%                     figure(100)
%     ax(1)= subplot(2,1,1)
%     imagesc(tout,fr_temp(good_freq),dyadic_spectrograms_temp(good_freq,:));%bsxfun(@rdivide,zsc_dy,sum(zsc_dy,1)))
%    % caxis([0 0.1])
%     ax(2)= subplot(2,1,2)
%     plot(0:1./300:(numel(ML_features.joint_angles.(angle_lists{kk}{ll}))-1)./300,ML_features.joint_angles.(angle_lists{kk}{ll})')
%     linkaxes(ax,'x')

%
% figure(99)
% plot(fr_temp(good_freq),bsxfun(@times,dyadic_spectrograms_temp(good_freq,1:50:500),1))
% %plot(dct(dyadic_spectrograms_temp(good_freq,1:50:1000)))
% %figure(99)
% figure(97)
% plot(fr_temp(good_freq),var(dyadic_spectrograms_temp(good_freq,:),[],2)./mean(dyadic_spectrograms_temp(good_freq,:),2))
%  %plot(fr_temp(good_freq),zsc_dy(good_freq,1:50:10000))
% %


%     [dyadic_spectrograms_temp,fr_temp,tout] = get_dyadic_spectrogram(ML_features.joint_angles.(angle_lists{kk}{ll})',opts);
%     zsc_dy = exp(dyadic_spectrograms_temp(:,:));
%     zsc_dy = bsxfun(@rdivide,bsxfun(@minus,zsc_dy,median(zsc_dy,2)),mad(zsc_dy,[],2));%,mad(zsc_dy,[],2));%zsc_dy
%     % zsc_dy = exp(dyadic_spectrograms_temp(:,:));
%     % zsc_dy(zsc_dy>1) = 1;
%     zsc_dy(zsc_dy>5) = 5;
%     zsc_dy(zsc_dy<-1) = -1;
% good_freq = find(fr_temp<30);


%% tsne garbage
% % % % % plot3(ML_features.spectrogram_pcs_hipR(:,1),ML_features.spectrogram_pcs_hipR(:,2),ML_features.spectrogram_pcs_hipR(:,3),'+')
% % % % %
% % % % %
%    subset = 1:100:100*2400;
%   mapped_trunk = tsne(cat(2,ML_features.spectrogram_pcs_wl_trunk_angle(subset,1:15)));
%     mapped_head = tsne(cat(2,ML_features.spectrogram_pcs_wl_trunk_angle(subset,1:15)));
%    mapped_jt_dynamics = tsne(cat(2,ML_features.spectrogram_pcs_wl_head_angle(subset,1:15),ML_features.spectrogram_pcs_wl_trunk_angle(subset,1:15)));
%       mapped_jt_dynamics_angles = tsne(cat(2,ML_features.spectrogram_pcs_wl_head_angle(subset,1:15),ML_features.spectrogram_pcs_wl_trunk_angle(subset,1:15)));
%
% % % % %
%    figure(22)
% plot(mapped(:,1),mapped(:,2),'+')

%  figure(102)
%  imagesc(ML_features.spectrogram_coeffs_head)
%


% get appearance feature
%ML_features.eigenpose_velocity
%ML_features.eigenpose_sd


%get_filtered_derivative(pose_score_whitened(:,1),numframes);
%get_window_sd(pose_score_whitened(:,1),numframes);

% get marker relative window features
%get_markersubgroup_velocity(mocapstruct.markers_preproc,[4:8],params);
%get_filtered_derivative(pose_score_whitened(:,1),numframes);
%get_window_sd(pose_score_whitened(:,1),numframes);

%     %  [dyadic_spectrograms,fr,~] = get_dyadic_spectrogram( agg_features_here',opts);
%     dyadic_spectrograms = [];
%     for jj = 1:size(agg_features_here,2)
% %                [~,fr_temp,~,dyadic_spectrograms_temp] = spectrogram(agg_features_here(:,jj),opts.clustering_window,...
% %           opts.clustering_overlap,1:30,opts.fps);
%
% %         params.tapers = [3 5];
% %         params.Fs = 300;
% %
% %         [S,t,f]=mtspecgramc( bsxfun(@minus,agg_features_here(:,jj),mean(agg_features_here(:,jj),1)),[1 1],params);
% %         dyadic_spectrograms_temp = (S(:,1:50).^2)';
% %
%
%
%
%        figure(102)
%        imagesc(1:size(dyadic_spectrograms_temp,2),fr_temp,(dyadic_spectrograms_temp))
%
%                 imagesc(log(dyadic_spectrograms_temp))
%
%       %  dyadic_spectrograms_temp = bsxfun(@rdivide,dyadic_spectrograms_temp,sum(dyadic_spectrograms_temp,2));
%         ML_features.mean_ja_spect{kk}(jj,:) = mean(log(dyadic_spectrograms_temp),2)';
%
%
%         dyadic_spectrograms = cat(1,dyadic_spectrograms,dyadic_spectrograms_temp);
%     end
%     dyadic_spectrograms = log( dyadic_spectrograms);

%          [dyadic_spectrograms_temp,tout,fr_temp]=mtspecgramc( bsxfun(@minus,agg_features_here(:,jj),mean(agg_features_here(:,jj),1)),[1 0.5],params);
%          dyadic_spectrograms_temp = dyadic_spectrograms_temp.^2;
%          dyadic_spectrograms_temp = log(dyadic_spectrograms_temp)';
%         dyadic_spectrograms_temp = dyadic_spectrograms_temp.^2;
%         % [~,fr_temp,~,dyadic_spectrograms_temp] = spectrogram(agg_features_here(:,jj),opts.clustering_window,...
%         %  opts.clustering_overlap,1:30,opts.fps);
%         dyadic_spectrograms = cat(1,dyadic_spectrograms,dyadic_spectrograms_temp);
%


% zsc_dy = (dyadic_spectrograms_temp(:,:));
% zsc_dy = bsxfun(@times,bsxfun(@rdivide,bsxfun(@minus,zsc_dy,median(zsc_dy,2)),mad(zsc_dy,[],2)),sqrt(zsc_dy));%,mad(zsc_dy,[],2));%zsc_dy
% zsc_dy = exp(dyadic_spectrograms_temp(:,:));
% zsc_dy(zsc_dy>1) = 1;
%zsc_dy = bsxfun(@rdivide,zsc_dy,sum(zsc_dy,1));
%     zsc_dy(zsc_dy>5) = 5;
%     zsc_dy(zsc_dy<-1) = -1;

%
% figure(99)
% plot(fr_temp(good_freq),bsxfun(@times,dyadic_spectrograms_temp(good_freq,1:50:1000),weighting_function))
% plot(dct(dyadic_spectrograms_temp(good_freq,1:50:1000)))
% figure(99)
%
%  plot(fr_temp(good_freq),zsc_dy(good_freq,1:50:10000))
%
%    figure(100)
%     ax(1)= subplot(2,1,1)
%     imagesc(tout,fr_temp(good_freq),bsxfun(@times,dyadic_spectrograms_temp(good_freq,:),weighting_function));%bsxfun(@rdivide,zsc_dy,sum(zsc_dy,1)))
%    caxis([0 0.1])
%     ax(2)= subplot(2,1,2)
%     plot(0:1./300:(size(ag

