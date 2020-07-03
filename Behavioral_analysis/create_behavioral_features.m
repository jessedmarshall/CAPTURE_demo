function ML_features = create_behavioral_features(mocapstruct,coeff_file,overwrite_coeff)

%% compute the joint angles
ML_features = compute_joint_angles_demo(mocapstruct);
%% compute the principal components of the joint angles and
ML_features = compute_appendage_pc_demos(mocapstruct,ML_features,coeff_file,overwrite_coeff);
%% compute the wavelet transform
ML_features = compute_wl_transform_features_demo(mocapstruct,ML_features,coeff_file,overwrite_coeff);

%% add window/vel/old features



%% code for visualization