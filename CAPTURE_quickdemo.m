%% Demo for CAPTURE Data
%[fList,pList] = matlab.codetools.requiredFilesAndProducts('CAPTURE_quickdemo.m');

% relies on following FEX/Open source contributions:
%Chronux
%Pca Randomized
%MTimesX
%Motionmapper
%Structvars

%load mocap file
datafile = ...
    load('C:\Users\Jesse Marshall\Documents\GitHub\Movement_analysis\Cortex_analysis\DemoRepo\Data\nolj_Recording_day8_caff1_nolj_imputed.mat');
mocapstruct = datafile;

%simple preprocess and turn into mocapstruct
%visualize the mocap data
animate_markers_nonaligned_fullmovie_demo(mocapstruct,1:10:10000);

% Create behavioral features
coefficient_file = 'demo_coefficients.mat';
overwrite_coefficient=0;
mocapstruct.modular_cluster_properties.clipped_index{8} = 1:size(mocapstruct.aligned_mean_position,1 );
MLmatobj = create_behavioral_features(mocapstruct,coefficient_file,overwrite_coefficient);

% perform a tsne embedding using a simple importance sampling
analysisparams.tsnegranularity = 50;
analysisstruct = compute_tsne_features(MLmatobj,mocapstruct,analysisparams);

zvals = tsne(analysisstruct.jt_features{8});
figure(1)
plot(zvals(:,1),zvals(:,2),'ob','MarkerFaceColor','b')

% or re-embed into Rat7M space

MLmatobj = create_extra_behavioral_features(mocapstruct,coefficient_file,overwrite_coefficient);

%visualize


%% run sequence and state analysis


%visualize