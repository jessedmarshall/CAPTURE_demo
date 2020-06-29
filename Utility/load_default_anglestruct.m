function anglestruct = load_default_anglestruct()


anglestruct.saggital_names = {'head_sagg','neck_sagg','spine_sagg'};
anglestruct.saggital_pairs =  {[2,3],[3,4],[4,5]}; %head, neck, spine angles , look in the z-y plane
%saggital_include = [1 1 1];

%transverse/overhead
anglestruct.transverse_names = {'head_trans','neck_trans','spine_trans','hipl_trans','hipr_trans',...
    'shouldl_trans','shouldr_trans','offset1_trans','offset2_trans'};
anglestruct.transverse_pairs =  {[2,3],[3,4],[4,5],[5,6],[5,7],[4 8 ], [4 9 ],[4,18],[4,19]}; %head, neck, spine angles , look in the z-y plane
%transverse_include = [1 1 1 1 1 0 0 1 1];

%coronal/along spine (front view)
anglestruct.coronal_names = {'head_coronal','hipl_coronal','hipr_coronal','shouldl_coronal','shouldr_coronal','offset1_coronal','offset2_coronal'};
anglestruct.coronal_pairs =  {[1,3],[5,6],[5,7],[4 8], [4 9],[4,18],[4,19]}; %head, neck, spine angles , look in the z-y plane
%coronal_include = [1 1 1 0 0 1 1];

%% angles to include for mai tsne
anglestruct.include_angles = {'head_sagg','neck_sagg','spine_sagg','head_trans','neck_trans',...
    'spine_trans','offset1_trans','offset2_trans','head_coronal','hipl_coronal','hipr_coronal',...
    'offset1_coronal','offset2_coronal','hipl_pitch','hipl_yaw','hipr_pitch','hipr_yaw'};% TO ADD: HIP YAW
%jointangle_struct.hipl_pitch,jointangle_struct.hipl_yaw

%alljt names
% allangles_names = {'lelbow_all','larm_all','relbow_all','rarm_all','lknee_all','lshin_all','rknee_all','rshin_all'};
% allangles_pairs =  {[8,10],[10,11],[9,12],[12,13],[6,14],[14,15],[7,16],[16,17]}; %head, neck, spine angles , look in the z-y plane
% all_include = zeros(1,numel(allangles_names));


%% get the various
anglestruct.segment_pairs = {{'HeadB','HeadL'},{'HeadF','HeadB'},{'HeadB','SpineF'},{'SpineF','SpineM'} ,...%1-4
    {'SpineL','SpineM'},{'SpineL','HipL'},{'SpineL','HipR'},... %5-7
    {'SpineF','ShoulderL'},{'SpineF','ShoulderR'},... %8,9
    {'ShoulderL','ElbowL'},{'ElbowL','ArmL'},{'ShoulderR','ElbowR'},{'ElbowR','ArmR'},...%10-13
    {'HipL','KneeL'},{'KneeL','ShinL'},{'HipR','KneeR'},{'KneeR','ShinR'},...%14-17
    {'SpineM','Offset1'},{'SpineM','Offset2'}}; %18,19


%hip angles to back of the spine
anglestruct.planar_trios{1}.plane = {{'SpineM','SpineL'},{'zvector','zvector'}};
anglestruct.planar_trios{1}.vector = {'SpineL','HipL'};
anglestruct.planar_trios{1}.name1 = 'hipl_pitch';
anglestruct.planar_trios{1}.name2 = 'hipl_yaw';
anglestruct.planar_trios{1}.namesuse = [1 2];

anglestruct.planar_trios{2}.plane = {{'SpineM','SpineL'},{'zvector','zvector'}};
anglestruct.planar_trios{2}.vector = {'SpineL','HipR'};
anglestruct.planar_trios{2}.name1 = 'hipr_pitch';
anglestruct.planar_trios{2}.name2 = 'hipr_yaw';
anglestruct.planar_trios{2}.namesuse = [1 2];
%knees in the plane defined by the hips and the back of the
%spine line
anglestruct.planar_trios{3}.plane = {{'SpineL','SpineM'},{'HipR','HipL'}};
anglestruct.planar_trios{3}.vector = {'HipL','KneeL'};
anglestruct.planar_trios{3}.name1 = 'kneel_yaw';
anglestruct.planar_trios{3}.name2 = 'kneel_pitch';
anglestruct.planar_trios{3}.namesuse = [1 2];

anglestruct.planar_trios{4}.plane = {{'SpineL','SpineM'},{'HipL','HipR'}};
anglestruct.planar_trios{4}.vector = {'HipR','KneeR'};
anglestruct.planar_trios{4}.name1 = 'kneer_yaw';
anglestruct.planar_trios{4}.name2 = 'kneer_pitch';
anglestruct.planar_trios{4}.namesuse = [1 2];


%shin line, relative to the plane of knee and hips, with axes
%chosen to get angles >0 to avoid ambiguity/0 crossing

anglestruct.planar_trios{5}.plane = {{'HipR','HipL'},{'HipL','KneeL'}};
anglestruct.planar_trios{5}.vector = {'KneeL','ShinL'};
anglestruct.planar_trios{5}.name1 = 'shinl_yaw';
anglestruct.planar_trios{5}.name2 = '~';
anglestruct.planar_trios{5}.namesuse = [1 ];

anglestruct.planar_trios{6}.plane = {{'HipL','HipR'},{'HipR','KneeR'}};
anglestruct.planar_trios{6}.vector = {'KneeR','ShinR'};
anglestruct.planar_trios{6}.name1 = 'shinr_yaw';
anglestruct.planar_trios{6}.name2 = '~';
anglestruct.planar_trios{6}.namesuse = [1 ];

anglestruct.planar_trios{7}.plane = {{'HipR','HipL'},{'HipL','KneeL'}};
anglestruct.planar_trios{7}.vector = {'KneeL','ShinL'};
anglestruct.planar_trios{7}.name1 = '~';
anglestruct.planar_trios{7}.name2 = 'shinl_pitch';
anglestruct.planar_trios{7}.namesuse = [ 2];

anglestruct.planar_trios{8}.plane = {{'HipL','HipR'},{'HipR','KneeR'}};
anglestruct.planar_trios{8}.vector = {'KneeR','ShinR'};
anglestruct.planar_trios{8}.name1 = '~';
anglestruct.planar_trios{8}.name2 = 'shinr_pitch';
anglestruct.planar_trios{8}.namesuse = [ 2];


% arm angles

anglestruct.planar_trios{9}.plane = {{'ShoulderR','ShoulderL'},{'ShoulderL','ElbowL'}};
anglestruct.planar_trios{9}.vector = {'ElbowL','ArmL'};
anglestruct.planar_trios{9}.name1 = 'arml_yaw';
anglestruct.planar_trios{9}.name2 = '~';
anglestruct.planar_trios{9}.namesuse = [1 ];

anglestruct.planar_trios{10}.plane = {{'ShoulderL','ShoulderR'},{'ShoulderR','ElbowR'}};
anglestruct.planar_trios{10}.vector = {'ElbowR','ArmR'};
anglestruct.planar_trios{10}.name1 = 'armr_yaw';
anglestruct.planar_trios{10}.name2 = '~';
anglestruct.planar_trios{10}.namesuse = [1 ];

anglestruct.planar_trios{11}.plane = {{'ShoulderL','ElbowL'},{'ShoulderR','ShoulderL'}};
anglestruct.planar_trios{11}.vector = {'ElbowL','ArmL'};
anglestruct.planar_trios{11}.name1 = '~';
anglestruct.planar_trios{11}.name2 = 'arml_pitch';
anglestruct.planar_trios{11}.namesuse = [ 2];

anglestruct.planar_trios{12}.plane = {{'ShoulderR','ElbowR'},{'ShoulderL','ShoulderR'}};
anglestruct.planar_trios{12}.vector = {'ElbowR','ArmR'};
anglestruct.planar_trios{12}.name1 = '~';
anglestruct.planar_trios{12}.name2 = 'armr_pitch';
anglestruct.planar_trios{12}.namesuse = [ 2];

%
%elbow angles

anglestruct.planar_trios{13}.plane = {{'SpineM','SpineF'},{'ShoulderR','ShoulderL'}};
anglestruct.planar_trios{13}.vector = {'ShoulderL','ElbowL'};
anglestruct.planar_trios{13}.name1 = 'elbowl_yaw';
anglestruct.planar_trios{13}.name2 = 'elbowl_pitch';
anglestruct.planar_trios{13}.namesuse = [1 2];

anglestruct.planar_trios{14}.plane = {{'SpineM','SpineF'},{'ShoulderL','ShoulderR'}};
anglestruct.planar_trios{14}.vector = {'ShoulderR','ElbowR'};
anglestruct.planar_trios{14}.name1 = 'elbowr_yaw';
anglestruct.planar_trios{14}.name2 = 'elbowr_pitch';
anglestruct.planar_trios{14}.namesuse = [1 2];
