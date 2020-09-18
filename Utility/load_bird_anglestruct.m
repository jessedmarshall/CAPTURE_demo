function anglestruct = load_bird_anglestruct()
    birdlinks = load('X:\Jesse\MotionAnalysisCaptures\DANNCE_animals\manuscript_formattedData\bird\bird18.mat');
         links =birdlinks.joints_idx;
birdlinks.joint_names{4} = 'SpineF';
birdlinks.joint_names{6} = 'SpineM';

%% use the segments in the skeleton
for lk = 1:size(links,1)
  anglestruct.segment_pairs{lk} =   {birdlinks.joint_names{links(lk,1)},birdlinks.joint_names{links(lk,2)}};
end



anglestruct.saggital_names = {'beak_sagg','head_sagg','spine_f_sagg','spine_r_sagg','chest_sagg','tail_sagg',...
    'lthigh_sagg','lankle_sagg','rthigh_sagg','rankle_sagg','lwing_sagg','rwing_sagg'};
anglestruct.saggital_pairs =  {[1,2],[2,3],[3,4],[4,6],[4,5],[6,7],...
    [6,12],[12,13],[6,17],[17,18],[10,11],[15,16]};

anglestruct.transverse_names = {'beak_tran','head_tran','spine_f_tran','spine_r_tran','chest_tran','tail_tran',...
  'leye_tran','lneck_tran','lwing_tran','reye_tran','rneck_tran','rwing_tran'};
anglestruct.transverse_pairs =  {[1,2],[2,3],[3,4],[4,6],[4,5],[6,7],...
[2,9],[9,10],[10,11],[2,14],[14,15],[15,16]};

anglestruct.coronal_names = {'beak_coronal','head_coronal','spine_f_coronal','spine_r_coronal','chest_coronal','tail_coronal',...
  'leye_coronal','lneck_coronal','lwing_coronal','reye_coronal','rneck_coronal','rwing_coronal'};
anglestruct.coronal_pairs =  {[1,2],[2,3],[3,4],[4,6],[4,5],[6,7],...
[2,9],[9,10],[10,11],[2,14],[14,15],[15,16]};

%head, neck, spine angles , look in the z-y plane
%saggital_include = [1 1 1];


% anglestruct.saggital_names = strcat(birdlinks.joint_names,'_sagg')';
%  %head, neck, spine angles , look in the z-y plane
% %saggital_include = [1 1 1];
% 
% %transverse/overhead
% anglestruct.transverse_names = strcat(birdlinks.joint_names,'_trans')';
% %anglestruct.transverse_pairs =  {links}; %head, neck, spine angles , look in the z-y plane
% %transverse_include = [1 1 1 1 1 0 0 1 1];
% 
% %coronal/along spine (front view)
% anglestruct.coronal_names = strcat(birdlinks.joint_names,'_coronal')';
%anglestruct.coronal_pairs =  {links}; %head, neck, spine angles , look in the z-y plane
%coronal_include = [1 1 1 0 0 1 1];

%% angles to include for mai tsne
anglestruct.include_angles = {anglestruct.saggital_names,...
    anglestruct.transverse_names,anglestruct.coronal_names};% TO ADD: HIP YAW

  %  anglestruct.saggital_pairs{lk} =  links(lk,:);
   %     anglestruct.transverse_pairs{lk} =  links(lk,:);
  %  anglestruct.coronal_pairs{lk} =  links(lk,:);


 anglestruct.planar_trios = [];
% birdlinks.joint_names(links)
% %% get the various
% anglestruct.segment_pairs = {{'EarR','EarL'},{'Snout','EarR'},{'EarR','SpineF'},{'SpineF','SpineM'} ,...%1-4
%     {'tail_base_','SpineM'},{'tail_base_','HindlimbL'},{'tail_base_','HindlimbR'},... %5-7
%     {'SpineF','ForelimbL'},{'SpineF','ForelimbR'},... %8,9
%     {'ForelimbL','ForepawL'},{'ForelimbR','ForepawR'},...%10-13
%     {'HindlimbL','HindpawL'},{'HindlimbR','HindpawR'}...%14-17
%     }; %18,19

% 
% %hip angles to back of the spine
% anglestruct.planar_trios{1}.plane = {{'SpineM','tail_base_'},{'zvector','zvector'}};
% anglestruct.planar_trios{1}.vector = {'tail_base_','HindlimbL'};
% anglestruct.planar_trios{1}.name1 = 'HindlimbL_pitch';
% anglestruct.planar_trios{1}.name2 = 'HindlimbL_yaw';
% anglestruct.planar_trios{1}.namesuse = [1 2];
% 
% anglestruct.planar_trios{2}.plane = {{'SpineM','tail_base_'},{'zvector','zvector'}};
% anglestruct.planar_trios{2}.vector = {'tail_base_','HindlimbR'};
% anglestruct.planar_trios{2}.name1 = 'HindlimbR_pitch';
% anglestruct.planar_trios{2}.name2 = 'HindlimbR_yaw';
% anglestruct.planar_trios{2}.namesuse = [1 2];
% %knees in the plane defined by the hips and the back of the
% %spine line
% anglestruct.planar_trios{3}.plane = {{'tail_base_','SpineM'},{'HindlimbR','HindlimbL'}};
% anglestruct.planar_trios{3}.vector = {'HindlimbL','HindpawL'};
% anglestruct.planar_trios{3}.name1 = 'HindpawL_yaw';
% anglestruct.planar_trios{3}.name2 = 'HindpawL_pitch';
% anglestruct.planar_trios{3}.namesuse = [1 2];
% 
% anglestruct.planar_trios{4}.plane = {{'tail_base_','SpineM'},{'HindlimbL','HindlimbR'}};
% anglestruct.planar_trios{4}.vector = {'HindlimbR','HindpawR'};
% anglestruct.planar_trios{4}.name1 = 'HindpawR_yaw';
% anglestruct.planar_trios{4}.name2 = 'HindpawR_pitch';
% anglestruct.planar_trios{4}.namesuse = [1 2];
% 

%shin line, relative to the plane of knee and hips, with axes
%chosen to get angles >0 to avoid ambiguity/0 crossing

% anglestruct.planar_trios{5}.plane = {{'HindlimbR','HindlimbL'},{'HindlimbL','HindpawL'}};
% anglestruct.planar_trios{5}.vector = {'HindpawL','ShinL'};
% anglestruct.planar_trios{5}.name1 = 'shinl_yaw';
% anglestruct.planar_trios{5}.name2 = '~';
% anglestruct.planar_trios{5}.namesuse = [1 ];
% 
% anglestruct.planar_trios{6}.plane = {{'HindlimbL','HindlimbR'},{'HindlimbR','HindpawR'}};
% anglestruct.planar_trios{6}.vector = {'HindpawR','ShinR'};
% anglestruct.planar_trios{6}.name1 = 'shinr_yaw';
% anglestruct.planar_trios{6}.name2 = '~';
% anglestruct.planar_trios{6}.namesuse = [1 ];
% 
% anglestruct.planar_trios{7}.plane = {{'HindlimbR','HindlimbL'},{'HindlimbL','HindpawL'}};
% anglestruct.planar_trios{7}.vector = {'HindpawL','ShinL'};
% anglestruct.planar_trios{7}.name1 = '~';
% anglestruct.planar_trios{7}.name2 = 'shinl_pitch';
% anglestruct.planar_trios{7}.namesuse = [ 2];
% 
% anglestruct.planar_trios{8}.plane = {{'HindlimbL','HindlimbR'},{'HindlimbR','HindpawR'}};
% anglestruct.planar_trios{8}.vector = {'HindpawR','ShinR'};
% anglestruct.planar_trios{8}.name1 = '~';
% anglestruct.planar_trios{8}.name2 = 'shinr_pitch';
% anglestruct.planar_trios{8}.namesuse = [ 2];
% 
% 
% % arm angles
% 
% anglestruct.planar_trios{9}.plane = {{'ForelimbR','ForelimbL'},{'ForelimbL','ForepawL'}};
% anglestruct.planar_trios{9}.vector = {'ForepawL','ArmL'};
% anglestruct.planar_trios{9}.name1 = 'arml_yaw';
% anglestruct.planar_trios{9}.name2 = '~';
% anglestruct.planar_trios{9}.namesuse = [1 ];
% 
% anglestruct.planar_trios{10}.plane = {{'ForelimbL','ForelimbR'},{'ForelimbR','ForepawR'}};
% anglestruct.planar_trios{10}.vector = {'ForepawR','ArmR'};
% anglestruct.planar_trios{10}.name1 = 'armr_yaw';
% anglestruct.planar_trios{10}.name2 = '~';
% anglestruct.planar_trios{10}.namesuse = [1 ];
% 
% anglestruct.planar_trios{11}.plane = {{'ForelimbL','ForepawL'},{'ForelimbR','ForelimbL'}};
% anglestruct.planar_trios{11}.vector = {'ForepawL','ArmL'};
% anglestruct.planar_trios{11}.name1 = '~';
% anglestruct.planar_trios{11}.name2 = 'arml_pitch';
% anglestruct.planar_trios{11}.namesuse = [ 2];
% 
% anglestruct.planar_trios{12}.plane = {{'ForelimbR','ForepawR'},{'ForelimbL','ForelimbR'}};
% anglestruct.planar_trios{12}.vector = {'ForepawR','ArmR'};
% anglestruct.planar_trios{12}.name1 = '~';
% anglestruct.planar_trios{12}.name2 = 'armr_pitch';
% anglestruct.planar_trios{12}.namesuse = [ 2];

%
%elbow angles 
% 
% anglestruct.planar_trios{5}.plane = {{'SpineM','SpineF'},{'ForelimbR','ForelimbL'}};
% anglestruct.planar_trios{5}.vector = {'ForelimbL','ForepawL'};
% anglestruct.planar_trios{5}.name1 = 'ForepawL_yaw';
% anglestruct.planar_trios{5}.name2 = 'ForepawL_pitch';
% anglestruct.planar_trios{5}.namesuse = [1 2];
% 
% anglestruct.planar_trios{6}.plane = {{'SpineM','SpineF'},{'ForelimbL','ForelimbR'}};
% anglestruct.planar_trios{6}.vector = {'ForelimbR','ForepawR'};
% anglestruct.planar_trios{6}.name1 = 'ForepawR_yaw';
% anglestruct.planar_trios{6}.name2 = 'ForepawR_pitch';
% anglestruct.planar_trios{6}.namesuse = [1 2];
