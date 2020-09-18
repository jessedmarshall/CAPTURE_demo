function ML_features = create_behavioral_features(mocapstruct,coeff_file,overwrite_coeff)
% Make features for tsne
% ---------------------------
% (C) Jesse D Marshall 2020
%     Harvard University 

%% compute the joint angles11


ML_features = compute_joint_angles_demo(mocapstruct);
%% compute the principal components of the joint angles and
ML_features = compute_appendage_pc_demos(mocapstruct,ML_features,coeff_file,overwrite_coeff);
%% compute the wavelet transform
tic
ML_features = compute_wl_transform_features_demo(mocapstruct,ML_features,coeff_file,overwrite_coeff);
toc
%% add window/vel/old features



%% code for visualization
