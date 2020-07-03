%% Demo for CAPTURE Data
% This repository contains contributions from the following FEX/Open source contributions which are included:
%Chronux
%Pca Randomized
%MTimesX
%othercolor
%Motionmapper
%Structvars
% ---------------------------
% (C) Jesse D Marshall 2020
%     Harvard University 





%load mocap file
datafile = ...
    load('C:\Users\Jesse Marshall\Documents\GitHub\Movement_analysis\Cortex_analysis\DemoRepo\Data\nolj_Recording_day8_caff1_nolj_imputed.mat');
mocapstruct = datafile;

%visualize the mocap data
animate_markers_nonaligned_fullmovie_demo(mocapstruct,1:10:10000);

% Create behavioral features
coefficient_file = 'demo_coefficients.mat';
overwrite_coefficient=0;
mocapstruct.modular_cluster_properties.clipped_index{8} = 1:size(mocapstruct.aligned_mean_position,1 );
MLmatobj = create_behavioral_features(mocapstruct,coefficient_file,overwrite_coefficient);

% perform a tsne embedding subselecting every 50 frames
analysisparams.tsnegranularity = 50;

%subselect a particular set of features
analysisstruct = compute_tsne_features(MLmatobj,mocapstruct,analysisparams);

%run tsne
zvals = tsne(analysisstruct.jt_features);
figure(1)
plot(zvals(:,1),zvals(:,2),'ob','MarkerFaceColor','b')
analysisstruct.zValues = zvals;

%% clustering parameters
analysisstruct.params.density_res = 1001; %resolution of the map
analysisstruct.params.density_width = 2; %density kernel in tsne space
analysisstruct.params.expansion_factor = 1.1; %add a little room to the map after kernel smoothing
analysisstruct.params.density_threshold = 1*10^(-5); %remove regions in plots with low density
analysisstruct.matchedconds = {[unique(analysisstruct.condition_inds)]}; %if running over multiple conditions
analysisstruct.conditions_to_run = [unique(analysisstruct.condition_inds)];
analysisstruct.tsnegranularity = analysisparams.tsnegranularity;
params.reorder=1;
analysisstruct = compute_analysis_clusters_demo(analysisstruct,params);

%% behavior plots and movies
analysisstruct.conditionnames = 'myrat';
analysisstruct.ratnames = 'myrat';
analysisstruct.filesizes = {540000};

h1=figure(609)
clf;
params.nameplot=0;
params.density_plot =0;
params.watershed = 1;
params.sorted = 1;
params.markersize = 1;
params.coarseboundary =0;
params.do_coarse = 0;
plot_clustercolored_tsne(analysisstruct,1,params.watershed,h1,params)
set(gcf,'Position',([100 100 1100 1100]))

animate_markers_nonaligned_fullmovie_demo(analysisstruct.mocapstruct_reduced_agg{1},...
    find(analysisstruct.annot_reordered{end}==100));

%% or use extnded set of 140 features
savefilename ='myextratsnefeature';
MLmatobj_extra = create_extra_behavioral_features(mocapstruct,'myrat',savefilename,overwrite_coefficient);
jt_features_extra = load_extra_tsne_features(mocapstruct,MLmatobj_extra,analysisparams);

% look at tsne of these added features
zvals_extra = tsne(jt_features_extra);
% or the combination
%zvals_extra = tsne(cat(2,analysisstruct.jt_features,jt_features_extra));
figure(2)
plot(zvals_extra(:,1),zvals_extra(:,2),'ob','MarkerFaceColor','b')
analysisstruct.zValues_extra = zvals;



%% run sequence and state analysis
params.do_show_pdistmatrix =1;
params.decimation_factor = 5; %downsample if needed to save on memory
params.doclustering = 1;

%clustering parameters
params.corr_threshold = 0.2;
params.clustercutoff = 0.65;
analysisstruct.plotdirectory = '';
%timescale to use, in seconds
params.timescales = [1./4 2]; 

analysisstruct.conditionnames = {'test'};
analysisstruct.ratname = {'myrat'};

hierarchystruct=   find_sequences_states_demo(analysisstruct,{'test'},params);

animate_markers_nonaligned_fullmovie_demo(analysisstruct.mocapstruct_reduced_agg{1},...
    find(hierarchystruct.clustered_behavior{1}==2));
%visualize


%internal: check dependencies
%[fList,pList] = matlab.codetools.requiredFilesAndProducts('CAPTURE_quickdemo.m');
