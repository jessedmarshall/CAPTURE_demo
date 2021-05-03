
basefolder = 'Y:\Jesse\Data\DeepLabCut_comparisons\dlc_megatrain_DATASET\';


datasetfile = load(strcat(basefolder,'dataset_new.mat'));

file_example = 50000;
image_file = strrep(strrep(datasetfile.dataset(file_example).image,'../',basefolder),'/',filesep);

%cameras 1-6 are CameraR, CameraL, CameraE, CameraU, CameraS, CameraU2.
filedays = {'_e0','_e1','_e2','_e3'};
fileday_ind = zeros(1,4);
for jj=1:4
fileday_ind(jj)=numel(strfind(image_file,filedays{jj}));
end


cameranames = {'CameraR','CameraL','CameraE','CameraU','CameraS','CameraU2'};
camera_ind = zeros(1,6);
for jj=1:6
camera_ind(jj)=numel(strfind(image_file,cameranames{jj}));
end

filetag = filedays{find(fileday_ind)};
threed_data = load(strcat(basefolder,'data_3d',filetag,'\',filetag(end),'_cam',num2str(find(camera_ind)),'_data.mat'));
calibration_data = load(strcat(basefolder,'calibration',filetag,...
    '\hires_cam',num2str(find(camera_ind)),'_params_rRDistort.mat'));

image_file_reduced = strrep(image_file,strcat('_',cameranames{find(camera_ind)},'.png'),'');
image_file_reduced = strrep(image_file_reduced,strcat(basefolder,'rat7M',...
    filedays{find(fileday_ind)},filesep,'sample'),'');
timeind = str2num(image_file_reduced(3:end));
%% find the correct sample
dataindex = find(threed_data.data_sampleID == timeind);

threed_sample = reshape(threed_data.data_3d(dataindex,:),3,20);
 twod_sample = reshape(threed_data.data_2d(dataindex,:),2,20)';


cameraParams = cameraParameters('IntrinsicMatrix',calibration_data.K,...
    'RadialDistortion',calibration_data.RDistort,...
    'TangentialDistortion',calibration_data.TDistort...
    ); 

imagePoints = worldToImage(cameraParams,...
    calibration_data.r,calibration_data.t,threed_sample','ApplyDistortion',true);

%% either run the
imagehere = imread(image_file);

figure(33)
clf;
imagesc(imagehere);
for  pt_ind = 1:20
    hold on
    plot(imagePoints(pt_ind,1)-21,imagePoints(pt_ind,2)-1,'ok','MarkerFaceColor','k')
    % slight rounding errors between different calc methods
     %   plot( datasetfile.dataset(file_example).joints{1}(pt_ind,2),...
      %      datasetfile.dataset(file_example).joints{1}(pt_ind,3),'or','MarkerFaceColor','r')
       %         plot( twod_sample(pt_ind,1)-21,...
       %     twod_sample(pt_ind,2)-1,'og','MarkerFaceColor','g')
end
hold off