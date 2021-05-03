function anglestruct = load_mouse_anglestruct()


anglestruct.saggital_names = {'head_sagg','neck_sagg','spine_sagg'};
anglestruct.saggital_pairs =  {[2,3],[3,4],[4,5]}; %head, neck, spine angles , look in the z-y plane
%saggital_include = [1 1 1];

%transverse/overhead
anglestruct.transverse_names = {'head_trans','neck_trans','spine_trans','HindlimbL_trans','HindlimbR_trans',...
    'shouldl_trans','shouldr_trans','offset1_trans','offset2_trans'};
anglestruct.transverse_pairs =  {[2,3],[3,4],[4,5],[5,6],[5,7],[4 8 ], [4 9 ]}; %head, neck, spine angles , look in the z-y plane
%transverse_include = [1 1 1 1 1 0 0 1 1];

%coronal/along spine (front view)
anglestruct.coronal_names = {'head_coronal','HindlimbL_coronal','HindlimbR_coronal','shouldl_coronal','shouldr_coronal','offset1_coronal','offset2_coronal'};
anglestruct.coronal_pairs =  {[1,3],[5,6],[5,7],[4 8], [4 9]}; %head, neck, spine angles , look in the z-y plane
%coronal_include = [1 1 1 0 0 1 1];

%% angles to include for mai tsne
anglestruct.include_angles = {'head_sagg','neck_sagg','spine_sagg','head_trans','neck_trans',...
    'spine_trans','offset1_trans','offset2_trans','head_coronal','HindlimbL_coronal','HindlimbR_coronal',...
    'offset1_coronal','offset2_coronal','HindlimbL_pitch','HindlimbL_yaw','HindlimbR_pitch','HindlimbR_yaw'};% TO ADD: HIP YAW
%jointangle_struct.HindlimbL_pitch,jointangle_struct.HindlimbL_yaw

%alljt names
% allangles_names = {'lelbow_all','larm_all','relbow_all','rarm_all','lknee_all','lshin_all','rknee_all','rshin_all'};
% allangles_pairs =  {[8,10],[10,11],[9,12],[12,13],[6,14],[14,15],[7,16],[16,17]}; %head, neck, spine angles , look in the z-y plane
% all_include = zeros(1,numel(allangles_names));


%% get the various
anglestruct.segment_pairs = {{'EarR','EarL'},{'Snout','EarR'},{'EarR','SpineF'},{'SpineF','SpineM'} ,...%1-4
    {'tail_base_','SpineM'},{'tail_base_','HindlimbL'},{'tail_base_','HindlimbR'},... %5-7
    {'SpineF','ForelimbL'},{'SpineF','ForelimbR'},... %8,9
    {'ForelimbL','ForepawL'},{'ForelimbR','ForepawR'},...%10-13
    {'HindlimbL','HindpawL'},{'HindlimbR','HindpawR'}...%14-17
    }; %18,19


%hip angles to back of the spine
anglestruct.planar_trios{1}.plane = {{'SpineM','tail_base_'},{'zvector','zvector'}};
anglestruct.planar_trios{1}.vector = {'tail_base_','HindlimbL'};
anglestruct.planar_trios{1}.name1 = 'HindlimbL_pitch';
anglestruct.planar_trios{1}.name2 = 'HindlimbL_yaw';
anglestruct.planar_trios{1}.namesuse = [1 2];

anglestruct.planar_trios{2}.plane = {{'SpineM','tail_base_'},{'zvector','zvector'}};
anglestruct.planar_trios{2}.vector = {'tail_base_','HindlimbR'};
anglestruct.planar_trios{2}.name1 = 'HindlimbR_pitch';
anglestruct.planar_trios{2}.name2 = 'HindlimbR_yaw';
anglestruct.planar_trios{2}.namesuse = [1 2];
%knees in the plane defined by the hips and the back of the
%spine line
anglestruct.planar_trios{3}.plane = {{'tail_base_','SpineM'},{'HindlimbR','HindlimbL'}};
anglestruct.planar_trios{3}.vector = {'HindlimbL','HindpawL'};
anglestruct.planar_trios{3}.name1 = 'HindpawL_yaw';
anglestruct.planar_trios{3}.name2 = 'HindpawL_pitch';
anglestruct.planar_trios{3}.namesuse = [1 2];

anglestruct.planar_trios{4}.plane = {{'tail_base_','SpineM'},{'HindlimbL','HindlimbR'}};
anglestruct.planar_trios{4}.vector = {'HindlimbR','HindpawR'};
anglestruct.planar_trios{4}.name1 = 'HindpawR_yaw';
anglestruct.planar_trios{4}.name2 = 'HindpawR_pitch';
anglestruct.planar_trios{4}.namesuse = [1 2];


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

anglestruct.planar_trios{5}.plane = {{'SpineM','SpineF'},{'ForelimbR','ForelimbL'}};
anglestruct.planar_trios{5}.vector = {'ForelimbL','ForepawL'};
anglestruct.planar_trios{5}.name1 = 'ForepawL_yaw';
anglestruct.planar_trios{5}.name2 = 'ForepawL_pitch';
anglestruct.planar_trios{5}.namesuse = [1 2];

anglestruct.planar_trios{6}.plane = {{'SpineM','SpineF'},{'ForelimbL','ForelimbR'}};
anglestruct.planar_trios{6}.vector = {'ForelimbR','ForepawR'};
anglestruct.planar_trios{6}.name1 = 'ForepawR_yaw';
anglestruct.planar_trios{6}.name2 = 'ForepawR_pitch';
anglestruct.planar_trios{6}.namesuse = [1 2];
