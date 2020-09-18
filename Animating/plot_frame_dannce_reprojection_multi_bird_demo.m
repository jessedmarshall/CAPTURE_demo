% %test reprojections wiht dannce
% timbase = '/home/twd/Dropbox/';
% basefolder = 'E:\Dropbox (Olveczky)\JesseMarshall\';
% timfilesep = '/';
%
% predictionstruct = load('X:\Jesse\MotionAnalysisCaptures\DANNCE_animals\test\predictions.mat');

function M=plot_frame_dannce_reprojection_multi_bird_demo(predictionstruct,frames_to_plot_vid_in,basefolder,camerause)


if nargin<3
    basefolder = 'E:\Dropbox (Olveczky)\JesseMarshall\';
    camerause=1;
end



if ~isfield(predictionstruct,'camnames_agg')
    predictionstruct.camnames_agg{1} = predictionstruct.cameras;
    predictionstruct.individual_index_vid = ones(1,size(predictionstruct.predictions.SpineM,1));
    predictionstruct.individual_index = ones(1,size(predictionstruct.predictions.SpineM,1));
end
frames_to_plot_vid_in_orig = frames_to_plot_vid_in;
unique_mice = unique(predictionstruct.individual_index_vid(frames_to_plot_vid_in));
for nn=unique_mice
    condframes = find(predictionstruct.individual_index_vid==nn);
    condframes_pts = find(predictionstruct.individual_index==nn);
    
    frames_to_plot_vid_in = find(ismember(condframes,frames_to_plot_vid_in_orig));
    
    if numel( frames_to_plot_vid_in)<numel(frames_to_plot_vid_in_orig)
        frames_to_plot_vid_in = cat(2,frames_to_plot_vid_in,frames_to_plot_vid_in(end)*ones(1,numel(frames_to_plot_vid_in_orig)-numel( frames_to_plot_vid_in)));
    end
    %
    % %test reprojections wiht dannce
    % timbase = '/home/twd/Dropbox/';
    % timfilesep = '/';
    %
    % videodirectory = strrep(strrep(predictionstruct.fullpath,timbase,...
    %     basefolder ),timfilesep,filesep);
    
    cameranames = fieldnames(predictionstruct.camnames_agg{nn});
    cameramatrix = cell(1,numel(cameranames));
    cameradirectory = cell(1,numel(cameranames));
    
    for lk=1:numel(cameranames)
        K=predictionstruct.camnames_agg{nn}.(cameranames{lk}).IntrinsicMatrix;
        RDistort=predictionstruct.camnames_agg{nn}.(cameranames{lk}).RadialDistortion;
        TDistort=double(predictionstruct.camnames_agg{nn}.(cameranames{lk}).TangentialDistortion);
        cameramatrix{lk} = cameraParameters('IntrinsicMatrix',K,'RadialDistortion',RDistort,'TangentialDistortion',TDistort);
        %  cameradirectory_base=   clean_tim_cameraname(predictionstruct.camnames_agg{nn}.(cameranames{lk}).video_directory);
        
        %% YOu will need to edit with your appropriate path
        cameradirectory_base = strrep(predictionstruct.camnames_agg{nn}.(cameranames{lk}).video_directory,...
            '/home/twd/Dropbox/autoencoder/selmaan/videos/',...
            'X:\Jesse\MotionAnalysisCaptures\DANNCE_animals\manuscript_formattedData\bird\videos\');
        
        
        
        %     cameradirectory_base =  ...
        %         strrep(strrep(predictionstruct.cameras.(cameranames{lk}).video_directory,timbase,...
        %         basefolder ),timfilesep,filesep);
        cameradirectory_cam = cameradirectory_base;
        %% flag to search if there is an xtra subdirectory
        dirnames = dir(cameradirectory_cam);
        %  if ~numel(strfind(cameradirectory_cam,'Camera'))
        if sum(cat(1,dirnames.isdir))>2 %% IF IN A SUBDIRECTORY
            %   cameradirectory_cam = strcat(cameradirectory_base,filesep,cameranames{lk});
            cameradirectory{lk} = strcat(cameradirectory_cam,filesep,dirnames(3).name);
        else
            cameradirectory{lk} =cameradirectory_cam;
        end
    end
    
    
    
    framename = strcat(cameradirectory{camerause},filesep,'0.avi');
    videoobj = VideoReader(framename);
    frameTimes = 1:max(predictionstruct.camnames_agg{nn}.(cameranames{camerause}).frame);
    frameTimes = double(double(frameTimes)./60);
    
    for frame_to_plot_vid = frames_to_plot_vid_in
        
        frame_to_plot_vid_new= double(predictionstruct.camnames_agg{nn}.(cameranames{camerause}).frame(frame_to_plot_vid));
        frame_to_plot_vid_true = frame_to_plot_vid_new;
        %get frames
        %
        % filebreak = 3500;
        % if numel(strfind(cameradirectory{camerause},'kyle'))
        %     predictionstruct.sample_factor=3;
        % filebreak = 3000;
        % elseif numel(strfind(cameradirectory{camerause},'marmoset'))
        % filebreak = 10000000;
        % end
        %
        %
        % fname = num2str(filebreak*floor(frame_to_plot_vid_new./filebreak));
        % frame_to_plot_vid_true = mod(frame_to_plot_vid_new,filebreak)+1;
        % framename = strcat(cameradirectory{camerause},filesep,fname,'.mp4');
        %
        %
        % if numel(strfind(cameradirectory{camerause},'marmoset'))
        %     framename = strrep(cameradirectory{camerause},'.avi','.mp4');
        %     frame_to_plot_vid_true = frame_to_plot_vid_new+1;
        %     framename = strrep(framename,'0.mp4','0_2.avi');
        %
        % end
        % frame_to_plot_points_mocap = frame_to_plot_vid*predictionstruct.sample_factor+predictionstruct.shift; %for points
        
        frame_to_plot_points = frame_to_plot_vid*predictionstruct.sample_factor+predictionstruct.shift; %for points
        frame_to_plot_points =condframes_pts(frame_to_plot_points);
        
        frame_to_plot_points(frame_to_plot_points<1) = 1;
        
        %get video -- don't open many video reader objects, which will eventually
        %crash matlab
        % if frame_to_plot_vid == frames_to_plot_vid_in(1)
        % framename_old = framename;
        % end
        %
        % if frame_to_plot_vid == frames_to_plot_vid_in(1) || ~strcmp(framename,framename_old)
        % videoobj = VideoReader(framename);
        % videoframe = read(videoobj,frame_to_plot_vid_true);%readFrame(videoobj);
        % %videoframe_undistorted = undistortImage(videoframe,cameramatrix{camerause});
        % videoframe_undistorted = videoframe;
        % framename_old = framename;
        % else
        %     videoframe = read(videoobj,frame_to_plot_vid_true);%readFrame(videoobj);
        %
        fprintf('frame time %f \n',frameTimes(frame_to_plot_vid_true))
        tic
        videoobj.CurrentTime= frameTimes(frame_to_plot_vid_true);
        toc
        videoframe = videoobj.readFrame;
        videoframe_undistorted = videoframe;
        %end
        
        
        %reproject
        points_to_plot = zeros(0,3);
        points_to_plot = [];
        markernames = fieldnames(predictionstruct.predictions);
        %markernames = markernames(1:(rend-1));
        for ll =1:numel(markernames)
            points_to_plot = cat(1,points_to_plot,predictionstruct.markers_preproc.(markernames{ll})(frame_to_plot_points,:));
        end
        points_to_plot = points_to_plot/525;
        
        imagePoints = worldToImage(cameramatrix{camerause},predictionstruct.camnames_agg{nn}.(cameranames{camerause}).rotationMatrix,...
            predictionstruct.camnames_agg{nn}.(cameranames{camerause}).translationVector,points_to_plot,'ApplyDistortion',true);
        
        %% look at reprojection in image from estimated world pose
        %figure(10+camerause)r
        set(gcf, 'Renderer', 'painters');
        imagesc(videoframe_undistorted)
        %clear videoobj
        
        hold on
        markercolors = predictionstruct.markercolor;
        numpoints =numel(markernames);
        
        for lk=1:(numpoints)% no tail tip for timesize(imagePoints,1)
            plot(imagePoints(lk,1),imagePoints(lk,2),'o','Color',markercolors{lk},'MarkerSize',4,'MarkerFaceColor',markercolors{lk})
        end
        hold off
        axis off
        
        
        %% add links
        hold on
        nmarkers = size(imagePoints,1);
        for mm = numel(predictionstruct.links):-1:1
            if numel(   predictionstruct.links{mm})
                if   predictionstruct.links{mm}(1)<=nmarkers && predictionstruct.links{mm}(2)<=nmarkers
                    if (ismember(predictionstruct.links{mm}(1),1:numel(predictionstruct.markernames)) && ismember(predictionstruct.links{mm}(2),1:numel(predictionstruct.markernames)))
                        % if (marker_plot(predictionstruct.links{mm}(1)) == 1 && marker_plot(predictionstruct.links{mm}(2)) == 1)
                        xx = [imagePoints(predictionstruct.links{mm}(1),1) ...
                            imagePoints(predictionstruct.links{mm}(2),1) ];
                        yy = [imagePoints(predictionstruct.links{mm}(1),2) ...
                            imagePoints(predictionstruct.links{mm}(2),2) ];
                        
                        
                        line(xx,yy,'Color',predictionstruct.markercolor{predictionstruct.links{mm}(1)},'LineWidth',2,'LineWidth',2,'Color',...
                            predictionstruct.markercolor{predictionstruct.links{mm}(1)})%,'MarkerFaceColor',predictionstruct.markercolor{predictionstruct.links{mm}(1)},'MarkerSize',6)
                        %end,'Marker','o'
                    end
                end
            end
        end
        
        hold off
        axis equal
        hold on
        
        for lk=1:(numpoints)% no tail tip for timesize(imagePoints,1)
            plot(imagePoints(lk,1),imagePoints(lk,2),'o','Color',markercolors{lk},'MarkerSize',5,'MarkerFaceColor',markercolors{lk})
        end
        hold off
        axis off
        
        
        
        %% for marm rotate camera angle
        if numel(strfind(cameradirectory{camerause},'marmoset'))
            if camerause == 2
                camroll(90);
            elseif camerause==3
                camroll(270);
            end
            
        end
        
        M(find(frame_to_plot_vid == frames_to_plot_vid_in)) =  getframe(gcf);
        clf;
    end
end
end