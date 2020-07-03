% preprocess data acquired using DANNCE
% ---------------------------
% (C) Jesse D Marshall 2020
%     Harvard University 

%[fList,pList] = matlab.codetools.requiredFilesAndProducts('preprocess_dannce.m');

% File Here to preprocess
datahere = load('C:\Users\Jesse Marshall\Documents\GitHub\Movement_analysis\Cortex_analysis\DemoRepo\Data\predictions.mat');

%do some surgery on names -- important for 
if isfield(datahere.predictions,'sampleID')
    datahere.sampleID = datahere.predictions.sampleID;
    datahere.predictions =rmfield(datahere.predictions,'sampleID');
end

% file specific changes in names
if isfield(datahere.predictions,'HeadBR')
    f = fieldnames(datahere.predictions);
    v = struct2cell(datahere.predictions);
    f{strmatch('HeadBR',f,'exact')} = 'HeadB';
    f{strmatch('HeadBL',f,'exact')} = 'HeadL';
    a = cell2struct(v,f);
    disp(a)
    datahere.predictions=a;
end

% parameters for preprocessing
preprocessing_parameters = struct();
preprocessing_parameters.median_filt_length = 5;
preprocessing_parameters.bad_frame_vel_thresh = 150; %this effectively turns the velocity criteria off
preprocessing_parameters.bad_frame_surround_flag = 0;
preprocessing_parameters.bad_frame_surround_number = 1;
preprocessing_parameters.interpolation_max_length = 5;
preprocessing_parameters.meanvelocity_lowpass = 60;
preprocessing_parameters.meanvelocity_lowpass = 60;
preprocessing_parameters.fastvelocity_threshold = 0.01;% 0.1;
preprocessing_parameters.moving_threshold = 0.00001;%0.015;
preprocessing_parameters.moving_framewindow = 600;

% the difference in framerate between the video and the canonical motion capture
% datasets
params.repfactor =10;

ratception_struct = preprocess_ratception_struct_demo(datahere,preprocessing_parameters,params);

%% load
[links,colors] = load_link_files('rats');
ratception_struct.links = links;
ratception_struct.markercolor = colors;
ratception_struct.markercolor = colors;

%% save
fileoutput = 'test_ratceptionstruct.mat';
save(fileoutput,'ratception_struct','-v7.3')

%%visualize
animate_markers_aligned_fullmovie_demo(ratception_struct,1:10:1000)
animate_markers_nonaligned_fullmovie_demo(ratception_struct,1:10:1000)

