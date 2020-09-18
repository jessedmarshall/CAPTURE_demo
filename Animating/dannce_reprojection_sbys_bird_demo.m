
function M=dannce_reprojection_sbys_bird_demo(predictionstruct,frames_to_plot_vid_in,basefolder,camerause,markercolors,numpoints)
if nargin<5
    markercolors = repmat({'r'},numel(predictionstruct.predictions));
end

if nargin<6
    numpoints =numel(markercolors);
end

if nargin<3
    basefolder = 'E:\Dropbox (Olveczky)\JesseMarshall\';
    camerause=1;
end


mathere = cell2mat( struct2cell( predictionstruct.markers_aligned_preproc));
maxper = prctile(mathere,97);
minper = prctile(mathere,3);


if ~isfield(predictionstruct,'camnames_agg')
    predictionstruct.camnames_agg{1} = predictionstruct.cameras;
        predictionstruct.individual_index_vid = ones(1,size(predictionstruct.predictions.SpineM,1));
        predictionstruct.individual_index = ones(1,size(predictionstruct.aligned_mean_position,1));
end
frames_to_plot_vid_in_orig = frames_to_plot_vid_in;
unique_mice = unique(predictionstruct.individual_index_vid(frames_to_plot_vid_in));
for nn=unique_mice
    condframes = find(predictionstruct.individual_index_vid==nn);
        condframes_pts = find(predictionstruct.individual_index==nn);

frames_to_plot_vid_in = find(ismember(condframes,frames_to_plot_vid_in_orig));
    %ismember(frames_to_plot_vid_in,predictionstruct.individual_index_vid)
cameranames = fieldnames(predictionstruct.camnames_agg{nn});
cameramatrix = cell(1,numel(cameranames));
cameradirectory = cell(1,numel(cameranames));

for lk=1:numel(cameranames)
    K=predictionstruct.camnames_agg{nn}.(cameranames{lk}).IntrinsicMatrix;
    RDistort=predictionstruct.camnames_agg{nn}.(cameranames{lk}).RadialDistortion;
    TDistort=double(predictionstruct.camnames_agg{nn}.(cameranames{lk}).TangentialDistortion);
    cameramatrix{lk} = cameraParameters('IntrinsicMatrix',K,'RadialDistortion',RDistort,'TangentialDistortion',TDistort);

    % cameradirectory_base = clean_tim_cameraname(predictionstruct.camnames_agg{nn}.(cameranames{lk}).video_directory);

    
    %% YOu will need to edit with your appropriate path
                camnamein = strrep(camnamein,'/home/twd/Dropbox/autoencoder/selmaan/videos/',...
                    'X:\Jesse\MotionAnalysisCaptures\DANNCE_animals\manuscript_formattedData\bird\videos\');
                
    cameradirectory_cam = cameradirectory_base;
        dirnames = dir(cameradirectory_cam);

     if ~numel(strfind(cameradirectory_cam,'dannce')) && numel(dirnames)>2%to exclude pup files for all but P7
 %   cameradirectory_cam = strcat(cameradirectory_base,filesep,cameranames{lk});
    cameradirectory{lk} = strcat(cameradirectory_cam,filesep,dirnames(1).name); %% change from 1 to 3 if an issue with pups
     elseif numel(strfind(cameradirectory_cam,'_pups')) && numel(dirnames)>2 %to exclude pup files?
  dirnames = dir(cameradirectory_cam);
    cameradirectory{lk} = strcat(cameradirectory_cam,filesep,dirnames(3).name);
     else
         %P7
        cameradirectory{lk} =cameradirectory_cam;
     end
    
     
    %dirnames = dir(cameradirectory_cam);
    %cameradirectory{lk} = strcat(cameradirectory_cam,filesep,dirnames(3).name);
end


    framename = strcat(cameradirectory{camerause},filesep,'0.avi');
videoobj = VideoReader(framename);
frameTimes = 1:max(predictionstruct.camnames_agg{nn}.(cameranames{camerause}).frame);
frameTimes = double(double(frameTimes)./60);

for frame_to_plot_vid = frames_to_plot_vid_in
    frame_to_plot_vid
  frame_to_plot_vid_new= double(predictionstruct.camnames_agg{nn}.(cameranames{camerause}).frame(frame_to_plot_vid));
      frame_to_plot_vid_true = frame_to_plot_vid_new;
% 
% %get frames
% if numel(strfind(cameradirectory{camerause},'marmoset'))
%     framename = strrep(cameradirectory{camerause},'.avi','.mp4');
%     frame_to_plot_vid_true = frame_to_plot_vid_new+1;
%     framename = strrep(framename,'0.mp4','0_2.avi');
% 
% elseif (strfind(cameradirectory{camerause},'kyle'))
%     fname = num2str(3000*floor(frame_to_plot_vid_new./3000));
% frame_to_plot_vid_true = mod(frame_to_plot_vid_new,3000)+1;
% framename = strcat(cameradirectory{camerause},filesep,fname,'.mp4');
% predictionstruct.sample_factor=3;
% elseif (strfind(cameradirectory{camerause},'bird'))
% else
% fname = num2str(3500*floor(frame_to_plot_vid_new./3500));
% frame_to_plot_vid_true = mod(frame_to_plot_vid_new,3500)+1;
% framename = strcat(cameradirectory{camerause},filesep,fname,'.mp4');
% end
       % frame_to_plot_points_mocap = frame_to_plot_vid*predictionstruct.sample_factor+predictionstruct.shift; %for points
%transform
    frame_to_plot_points = frame_to_plot_vid*predictionstruct.sample_factor+predictionstruct.shift; %for points
    
    frame_to_plot_points =condframes_pts(frame_to_plot_points);
    
frame_to_plot_points(frame_to_plot_points<1) = 1;


%% read by time for the bird

videoobj.CurrentTime= frameTimes(frame_to_plot_vid_true);
    videoframe = videoobj.readFrame;


%videoframe = read(videoobj,frame_to_plot_vid_true);%readFrame(videoobj);
% don't distort image
%videoframe_undistorted = undistortImage(videoframe,cameramatrix{camerause});
videoframe_undistorted = videoframe;

%reproject
points_to_plot = zeros(0,3);
points_to_plot = [];
markernames = fieldnames(predictionstruct.predictions);
markernames = markernames(1:(end));
for ll =1:numel(markernames)
     points_to_plot= cat(1,points_to_plot,predictionstruct.predictions.(markernames{ll})(frame_to_plot_points,:));
end
points_to_plot = points_to_plot/525;



imagePoints = worldToImage(cameramatrix{camerause},predictionstruct.camnames_agg{nn}.(cameranames{camerause}).rotationMatrix,...
    predictionstruct.camnames_agg{nn}.(cameranames{camerause}).translationVector,points_to_plot,'ApplyDistortion',true);
%% look at reprojection in image from estimated world pose
numpoints = numel(markernames);
msize = 2;
figure(10+camerause)
subplot(1,2,1)
imagesc(videoframe_undistorted)
hold on
for lk=1:numpoints%size(imagePoints,1)
    
plot(imagePoints(lk,1),imagePoints(lk,2),'o','Color',markercolors{lk},'MarkerSize',msize,'MarkerFaceColor',markercolors{lk})
end
hold off
axis off
axis equal

%% also plot links
  %% plot the links between markers
  hold on
  nmarkers = numpoints;
    for mm = numel(predictionstruct.links):-1:1
        if numel(   predictionstruct.links{mm})
        if   predictionstruct.links{mm}(1)<=nmarkers && predictionstruct.links{mm}(2)<=nmarkers
        if (ismember(predictionstruct.links{mm}(1),1:numel(predictionstruct.markernames)) && ismember(predictionstruct.links{mm}(2),1:numel(predictionstruct.markernames)))
           % if (marker_plot(predictionstruct.links{mm}(1)) == 1 && marker_plot(predictionstruct.links{mm}(2)) == 1)
                  xx = [imagePoints(predictionstruct.links{mm}(1),1) ...
                    imagePoints(predictionstruct.links{mm}(2),1) ];
                yy = [imagePoints(predictionstruct.links{mm}(1),2) ...
                    imagePoints(predictionstruct.links{mm}(2),2) ];
                
                   
                 line(xx,yy,'Color',predictionstruct.markercolor{predictionstruct.links{mm}(1)},'LineWidth',2,'Marker','o','LineWidth',2,'Color',...
                    predictionstruct.markercolor{predictionstruct.links{mm}(1)},'MarkerFaceColor',predictionstruct.markercolor{predictionstruct.links{mm}(1)},'MarkerSize',msize)
             %end
        end
        end
        end
    end

hold off






%% Plot wireframe
hh=subplot(1,2,2)
set(gcf,'Color','k')
predictionstruct.markernames = fieldnames(predictionstruct.predictions);
predictionstruct.markercolor = markercolors;
predictionstruct.mocapfiletimes = {[1]};
%predictionstruct.markers_preproc = predictionstruct.predictions;
%predictionstruct.markers_aligned_preproc = predictionstruct.predictions;
plot_frame_sub(predictionstruct,frame_to_plot_points,[],numpoints)
view([-8 13])

% zlim(2*[minper(:,3) maxper(:,3)])
% xlim(1.5*[minper(:,1) maxper(:,1)])
% ylim(1.5*[minper(:,2) maxper(:,2)])
    zlim([-110 150])
    xlim([-90 90])
    ylim([-90 90])
xlim
axis off
         M(find(frame_to_plot_vid == frames_to_plot_vid_in)) =  getframe(gcf);
end
end
end