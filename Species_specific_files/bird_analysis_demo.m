%% Demo for bird behavioral analysis
% ---------------------------
% (C) Jesse D Marshall 2020
%     Harvard University

%internal: check dependencies
[fList,pList] = matlab.codetools.requiredFilesAndProducts('bird_analysis_demo.m');


%% run the analysis
basedirectory = 'mydirectory';
%input predictions in DANNCE format
birdfilename = strcat(basedirectory,filesep,'predictions_new.mat');
%outputfile
birdfilename_out = strcat(basedirectory,filesep,'ratception_struct.mat');

input_params.SpineM_marker = 'centerBack';
input_params.SpineF_marker = 'backHead';
input_params.conversion_factor = 525; %mm/selman
input_params.repfactor = 5;

%% preprocess the data
ratception_struct = preprocess_dannce(birdfilename,birdfilename_out,'bird',input_params);

%% copy over camera information and metadata
predictionsfile = load(birdfilename);
if isfield(predictionsfile,'cameras')
    predictionsfieldnames = fieldnames(predictionsfile);
    for lk=1:numel(predictionsfieldnames)
        ratception_struct.(predictionsfieldnames{lk}) = predictionsfile.(predictionsfieldnames{lk});
    end
end
save(strcat(basedirectory,filesep,'ratception_struct.mat'),'-struct','ratception_struct','-v7.3')


%% visualize the wireframe in aligned coordinates
animate_markers_aligned_fullmovie_demo(ratception_struct,1:10:10000)
figure(370)
clf;

%% visualize in unaligned
% bird specific axes
axisparams.zlim = ([200 350]);
axisparams.xlim = ([-300 300]);
axisparams.zlim = ([-300 300]);

animate_markers_nonaligned_fullmovie_demo(ratception_struct,1:10:1000,[],axisparams)

% visualize reproject
%ratception_struct.cameras.lBack.video_directory
ratception_struct.predictions = ratception_struct.markers_preproc;
ratception_struct.sample_factor = 5;%bird is at 60 Hz
ratception_struct.shift = 0;

M=dannce_reprojection_sbys_bird_demo(ratception_struct,1:30:2000,...
    [],1,ratception_struct.markercolor,20)


%% do embedding
[analysisstruct,hierarchystruct] = CAPTURE_quickdemo(...
    strcat(basedirectory,filesep,'ratception_struct.mat'),...
    'bird','my_coefficients');

save(strcat(basedirectory,filesep,'myanalysisstruct.mat'),'-struct','analysisstruct',...
    '-v7.3')
save(strcat(basedirectory,filesep,'myhierarchystruct.mat'),'-struct','hierarchystruct',...
    '-v7.3')

%% plot the tsne
plotfolder = strcat(basedirectory,filesep,'plots\');
mkdir(plotfolder)

%% look at the tsne
h1=figure(608)
clf;
params.nameplot=1;
params.density_plot =0;
params.watershed = 1;
params.sorted = 1;
params.markersize = 0.2;
params.jitter = 0;
params.coarseboundary = 0;
analysisstruct.params.density_width=0.25;
analysisstruct.params.density_res=4001;
plot_clustercolored_tsne(analysisstruct,1,params.watershed,h1,params)
set(gcf,'renderer','painters')
colorbar off
axis equal
% set(gcf,'Position',([100 100 1100 1100]))
print('-dpng',strcat(plotfolder,'birdtsnev3.png'),'-r1200')
print('-depsc',strcat(plotfolder,'birdtsnev3.eps'),'-r1200')
%


%% go from the downsampled tsne to the full-res
maxframes = max(analysisstruct.frames_with_good_tracking{1})+(analysisstruct.tsnegranularity-1);
annot_frames_full = zeros(1,maxframes);
seq_frames_full = zeros(1,maxframes);
state_frames_full = zeros(1,maxframes);

frames_to_fill = unique(sort((reshape(bsxfun(@plus,analysisstruct.frames_with_good_tracking{1},[0:(analysisstruct.tsnegranularity-1)]),1,[]))));
annot_frames_full(frames_to_fill) = repelem(analysisstruct.annot_reordered{end}',analysisstruct.tsnegranularity,1);
seq_frames_full(frames_to_fill) = repelem(hierarchystruct.clustered_behavior{1},analysisstruct.tsnegranularity,1);
state_frames_full(frames_to_fill) = repelem(hierarchystruct.clustered_behavior{2},analysisstruct.tsnegranularity,1);




%% plot some sequences and states
animate_markers_aligned_fullmovie_demo(analysisstruct.mocapstruct_reduced_agg{1},...
    find(hierarchystruct.clustered_behavior{1}==4));

axisparams.zlim = ([100 350]);
axisparams.xlim = ([-400 400]);
axisparams.ylim = ([-400 400]);

seqname = {'sequences','state'};
for nn=1:2
    savefolder = strcat(basefolder,filesep,'output_videos',filesep,seqname{nn},filesep);
    mkdir(savefolder)
    for lk =1:max(hierarchystruct.clustered_behavior{nn})
        if nn==1
            indsplot = find(seq_frames_full==lk);
        else
            indsplot = find(state_frames_full==lk);
        end
        indsplot =indsplot(1:min(20000,numel(indsplot)));
        figure(370);
        clf;
        
        M=dannce_reprojection_sbys_bird_demo(ratception_struct,floor(indsplot(1:20:end)./5),...
            [],1,ratception_struct.markercolor,20)
        
        %% 2x RT
        vidfile_out = strcat(savefolder,filesep,'reprojection_v2_sequencenumber',num2str(lk),'.mp4');
        vwrite = VideoWriter(vidfile_out,'MPEG-4');
        vwrite.FrameRate = 30;
        
        open(vwrite)
        writeVideo(vwrite, M)
        close(vwrite)
    end
end






%% make individual videos per cluster
basefolder = 'X:\Jesse\MotionAnalysisCaptures\DANNCE_animals\manuscript_formattedData\bird\';

ratception_struct = load('X:\Jesse\MotionAnalysisCaptures\DANNCE_animals\manuscript_formattedData\bird\ratception_struct.mat');
ratception_struct.predictions = ratception_struct.markers_preproc;
analysisstruct = load(strcat(basefolder,filesep,'myanalysisstruct.mat'));

analysisstruct.mocapnames{1}{1} = matfile('X:\Jesse\MotionAnalysisCaptures\DANNCE_animals\manuscript_formattedData\bird\ratception_struct.mat');
clustuse=1;
[agg_mocap_structs,agg_snippetinds,agg_mocap_preproc] = collect_mocap_snippets_allzvals_demo(analysisstruct,1);
agg_mocap_structs_snippets = agg_mocap_structs;
agg_preproc = agg_mocap_preproc;
analysisstruct.agg_mocap_structs_snippets=agg_mocap_structs;
analysisstruct.agg_snippetinds=agg_snippetinds;
analysisstruct.agg_preproc=agg_mocap_preproc;


ratception_struct.sample_factor = 5;
ratception_struct.shift = 0;
ratception_struct.camuse = 1;

conditionhere=1;
mosaicfolder = strcat('X:\Jesse\ClusterVideos\DANNCE_videos_birdv2\');
mkdir(mosaicfolder)


%% do embedding
basefolder = 'X:\Jesse\MotionAnalysisCaptures\DANNCE_animals\manuscript_formattedData\bird';
ratception_struct.basefolder = basefolder;
%save_mosaicclusters(analysisstruct,1,mosaicfolder,2)

for chere =1%clustuse'
    figure(388)
    save_mosaicclusters_reproject_bird_demo(analysisstruct,chere,mosaicfolder,...
        1,[],ratception_struct)
end








