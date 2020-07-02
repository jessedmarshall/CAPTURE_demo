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

%subselect a particular set of features
analysisstruct = compute_tsne_features(MLmatobj,mocapstruct,analysisparams);

%run tsne
zvals = tsne(analysisstruct.jt_features);
figure(1)
plot(zvals(:,1),zvals(:,2),'ob','MarkerFaceColor','b')
analysisstruct.zValues = zvals;

%% cluster
analysisstruct.params.density_res = 1001;
analysisstruct.params.density_width = 2; %rat markerless is 1.5, mouse is 2 %kyle is    2
analysisstruct.params.expansion_factor = 1.1;
analysisstruct.params.density_threshold = 1*10^(-5);
analysisstruct.matchedconds = {[unique(analysisstruct.condition_inds)]};
analysisstruct.conditions_to_run = [unique(analysisstruct.condition_inds)];
analysisstruct.coarse_annotation_mat = [];
analysisstruct.tsnegranularity = analysisparams.tsnegranularity;
params.reorder=1;
analysisstruct = compute_analysis_clusters(analysisstruct,params);
%save analysisstruct


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

% or re-embed into Rat7M space
MLmatobj = create_extra_behavioral_features(mocapstruct,coefficient_file,overwrite_coefficient);


%% run sequence and state analysis
params.do_show_pdistmatrix =1;
params.decimation_factor = 5;
params.doclustering = 1;
%clustering parameters
params.corr_threshold = 0.2;
params.clustercutoff = 0.65;
analysisstruct.plotdirectory = '';
params.timescales = [1./4 2]; %timescale to use, in s

analysisstruct.conditionnames = {'test'};
analysisstruct.ratname = {'myrat'};

hierarchystruct=    find_sequences_states_demo(analysisstruct,{'test'},params);

animate_markers_nonaligned_fullmovie_demo(analysisstruct.mocapstruct_reduced_agg{1},...
    find(hierarchystruct.clustered_behavior{1}==2));
%visualize