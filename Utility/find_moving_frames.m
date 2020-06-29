function [moveframes,restframes,filteredvelocity,movementparams] = find_moving_frames(markers_preproc,preprocessing_parameters,showplots,velonly)
params.fps = 300;
params.difforder = 10;
params.medfiltorder = 3;
params.gaussorder = 2.5;
  params.difforder_movav = 300;
  showplots = 0;
  velonly = 0;
  
  fnames = fieldnames(markers_preproc);
  for ll = 1:numel(fnames)
    velocity_temp = sum(diff(markers_preproc.(fnames{ll}),1,1).^2,2);
    badfr = bsxfun(@plus,find(velocity_temp>15),(-1:1));
    badfr(badfr<1) = 1;
    badfr(badfr>size( markers_preproc.(fnames{ll}),1)) = size( markers_preproc.(fnames{ll}),1);
    markers_preproc.(fnames{ll})(badfr,:) = nan;
  end
  
  %get subgp velocity. should ignore nans
[av_vel_head,vel_comps,av_std,std_comps,av_accel,accel_comps,av_std_accel,std_comps_accel] = ...
    get_markersubgroup_velocity(markers_preproc,[1,2,3],params);
av_vel_head(isnan(av_vel_head)) = 0;

[av_vel_trunk,vel_comps,av_std,std_comps,av_accel,accel_comps,av_std_accel,std_comps_accel] = ...
    get_markersubgroup_velocity(markers_preproc,[4,5,6,7,8],params);
av_vel_trunk(isnan(av_vel_trunk)) = 0;

[av_vel_hips,vel_comps,av_std,std_comps,av_accel,accel_comps,av_std_accel,std_comps_accel] = ...
    get_markersubgroup_velocity(markers_preproc,[9,10],params);
av_vel_hips(isnan(av_vel_hips)) = 0;


sum_vel = sqrt(av_vel_head.^2+av_vel_trunk.^2+av_vel_hips.^2)./3;



%dHipass = designfilt('highpassiir', 'FilterOrder', 3, 'HalfPowerFrequency', 0.3/(params.fps/2), ...
 %   'DesignMethod', 'butter');
            
       
           % [f1,f2] = tf(dHipass);
          
                filteredvelocity =  sum_vel;%filtfilt(f1,f2, sum_vel);
                
if (nargin>1 && showplots==1)
figure(335)
plot(av_vel_head)
hold on
plot(av_vel_trunk,'r')
plot(av_vel_hips,'g')
plot(filteredvelocity,'k')
set(gca,'YScale','log')
ylim([0.01 100])

figure(334)
xspace = logspace(-3,2,100);
[n,x] = hist(filteredvelocity,xspace);
loglog(xspace,n)
end
if nargin <3 || velonly == 0
%thresh = 0.1;
thresh = preprocessing_parameters.fastvelocity_threshold;
fastframes = find(filteredvelocity>thresh);
nonfastframes = setxor(1:numel(filteredvelocity),fastframes);


movementparams = params;
movementparams.movethresh = thresh;

total_moveframes_fast = zeros(1,numel(filteredvelocity));
total_moveframes_fast(rectify_inds(unique(reshape(bsxfun(@plus,fastframes',(-preprocessing_parameters.moving_framewindow:...
    preprocessing_parameters.moving_framewindow)'),[],1)),numel(total_moveframes_fast))) = 1;
moveframes = find(total_moveframes_fast);
restframes = find(total_moveframes_fast==0);
else
    moveframes = [];
    restframes = [];
    movementparams =[];
end
%otherframes = setxor(mocapstruct_all.move_frames,intersect(mocapstruct_all.move_frames,total_moveframes_fast));
end

% 
% animate_markers_aligned_fullmovie(mocapstruct_all,indsfast(1:1000:end)')
% animate_markers_aligned_fullmovie(mocapstruct_all,indsrest(1:1000:end)')
% 
% animate_markers_aligned_fullmovie(mocapstruct_all,otherframes(1:300:end))
%animate_markers_aligned_fullmovie(mocapstruct_all,mocapstruct_all.rest_frames(1:1000:end))


