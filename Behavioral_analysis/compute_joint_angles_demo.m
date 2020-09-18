function ML_features = compute_joint_angles_demo(mocapstruct)
ML_features = struct();


%% get joint angle features
%saggital/cross section (ie side view)
fprintf('computing joing angles \n');

%% these are defined wrt the absolute coordinate system of the arena/animal
%saggital/pitch is angle in the plane of the long axis of the
%animal
%transverse/yaw is angle in the plane
%coronal is facing down the barrel of the animal
saggital_inds = [2,3];
coronal_inds = [1,3];
transverse_inds = [1,2];
allangles_inds = [1,2,3]; %use on knees and arms
%transverse_pairs

  % anglestruct = load_mouse_anglestruct() ;
  % anglestruct = load_mouse_kyle_anglestruct() ;
       anglestruct = load_default_anglestruct() ;
       anglestruct = load_bird_anglestruct() ;

segment_pairs=anglestruct.segment_pairs;
coronal_pairs=anglestruct.coronal_pairs;
saggital_pairs=anglestruct.saggital_pairs;
transverse_pairs=anglestruct.transverse_pairs;
planar_trios = anglestruct.planar_trios;
%% load the variables in from the structfield
   %anglestruct_fieldnames = fieldnames(anglestruct);
  assigns= structvars(anglestruct);
  for kk=1:size(assigns,1)
eval(  assigns(kk,:))
  end

jointangle_struct = struct();
%mean_seglengths = struct();
all_seglengths = cell(1,numel(segment_pairs));
all_segments = cell(1,numel(segment_pairs));
transverse_seglengths = cell(1,numel(segment_pairs));

ML_features.segment_pairs = segment_pairs;
ML_features.include_angles = include_angles;

for ll = 1:numel(saggital_pairs)
    vec1 =  mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(1)}{1})-...
        mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(1)}{2});
    
    vec2 =  mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(2)}{1})-...
        mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(2)}{2});

    jointangle_struct.(saggital_names{ll}) = ...
        acosd(dot(vec1(:,saggital_inds)',vec2(:,saggital_inds)')...
        ./(sqrt(sum(vec1(:,saggital_inds).^2,2)).*sqrt(sum(vec2(:,saggital_inds).^2,2)))')';

end

for ll = 1:numel(coronal_pairs)
    vec1 =  mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(1)}{1})-...
        mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(1)}{2});
    
    vec2 =  mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(2)}{1})-...
        mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(2)}{2});

    
    jointangle_struct.(coronal_names{ll}) = ...
        acosd(dot(vec1(:,coronal_inds)',vec2(:,coronal_inds)')...
        ./(sqrt(sum(vec1(:,coronal_inds).^2,2)).*sqrt(sum(vec2(:,coronal_inds).^2,2)))')'; 
end


%% transverse
for ll = 1:numel(transverse_pairs)
    vec1 =  mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(1)}{1})-...
        mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(1)}{2});
    
    vec2 =  mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(2)}{1})-...
        mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(2)}{2});
  
    jointangle_struct.(transverse_names{ll}) = ...
        acosd(dot(vec1(:,transverse_inds)',vec2(:,transverse_inds)')...
        ./(sqrt(sum(vec1(:,transverse_inds).^2,2)).*sqrt(sum(vec2(:,transverse_inds).^2,2)))')';

    
end

for ll = 1:numel(transverse_pairs)
     vec1 =  mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(1)}{1})-...
        mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(1)}{2});

    transverse_seglengths{ll} = (sqrt(sum(vec1(:,transverse_inds).^2,2)));
end


for kk=1:numel(planar_trios)
    if numel(planar_trios{kk}.namesuse)==2
[jointangle_struct.(planar_trios{kk}.name1),jointangle_struct.(planar_trios{kk}.name2)] =  get_planar_jointangles(mocapstruct,...
    planar_trios{kk}.plane, planar_trios{kk}.vector);
    elseif find(planar_trios{kk}.namesuse==1)
        [jointangle_struct.(planar_trios{kk}.name1),~] =  get_planar_jointangles(mocapstruct,...
    planar_trios{kk}.plane, planar_trios{kk}.vector);
    elseif find(planar_trios{kk}.namesuse==2)
        [~,jointangle_struct.(planar_trios{kk}.name2)] =  get_planar_jointangles(mocapstruct,...
    planar_trios{kk}.plane, planar_trios{kk}.vector);
    end
end




%% get the segment lengths
fprintf('getting segment lengths \n')
for ll = 1:numel(segment_pairs)
    vec1 =  mocapstruct.markers_aligned_preproc.(segment_pairs{ll}{1})-...
        mocapstruct.markers_aligned_preproc.(segment_pairs{ll}{2});
    
    all_seglengths{ll} = (sqrt(sum(vec1(:,allangles_inds).^2,2)));
    all_segments{ll} = vec1;
    % end
end


%total_include = cat(2,saggital_include,coronal_include,all_include,transverse_include);
%fn_ja_all = fieldnames(jointangle_struct);

%% get the elements that aren't directly organized in relation to body axes

%total_include = find(ismember(fieldnames(jointangle_struct),include_angles));

jointangle_struct = structfun(@(x) real(x),jointangle_struct, 'UniformOutput', false);
fname = fieldnames(jointangle_struct);
for lk=1:numel(fname)
    jointangle_struct.(fname{lk})(find(isnan( jointangle_struct.(fname{lk})))) = 0;
end
%mean_seglengths = structfun(@(x) real(x),mean_seglengths, 'UniformOutput', false);

ML_features.joint_angles_mean = real(structfun(@nanmean,(jointangle_struct)));
ML_features.jointangle_struct = jointangle_struct;
ML_features.all_seglengths = all_seglengths;
ML_features.all_segments = all_segments; 
ML_features.transverse_seglengths = transverse_seglengths; 


end

