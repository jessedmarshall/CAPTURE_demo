function [av_vel,vel_comps,av_std,std_comps,av_accel,accel_comps,av_std_accel,std_comps_accel] =...
    get_markersubgroup_velocity(markermocap,submarkers,params)
%shoould be sum not mean of components!!
fieldnames_here = fieldnames(markermocap);

num_components = size(markermocap.(fieldnames_here{(1)}),2);
%velocity
av_vel = zeros(size(markermocap.(fieldnames_here{(1)}),1),num_components);
vel_comps = zeros(size(markermocap.(fieldnames_here{(1)}),1),3);

%standard deviation
av_std = zeros(size(markermocap.(fieldnames_here{(1)}),1),num_components);
std_comps = zeros(size(markermocap.(fieldnames_here{(1)}),1),3);

%accel and accel std
av_accel = zeros(size(markermocap.(fieldnames_here{(1)}),1),num_components);
accel_comps = zeros(size(markermocap.(fieldnames_here{(1)}),1),3);

av_std_accel = zeros(size(markermocap.(fieldnames_here{(1)}),1),num_components);
std_comps_accel = zeros(size(markermocap.(fieldnames_here{(1)}),1),3);


for jj = 1:size(markermocap.(fieldnames_here{(1)}),2)
    vel_temp = zeros(size(markermocap.(fieldnames_here{(1)}),1),numel(submarkers));
        std_temp = zeros(size(markermocap.(fieldnames_here{(1)}),1),numel(submarkers));
         accel_temp = zeros(size(markermocap.(fieldnames_here{(1)}),1),numel(submarkers));
        accel_std_temp = zeros(size(markermocap.(fieldnames_here{(1)}),1),numel(submarkers));

    for mm = submarkers
        tracetemp = markermocap.(fieldnames_here{(mm)})(:,jj);
        tracetemp(tracetemp==0) = nan;
        [veloc,stand,accelout,accelstd] = get_vector_velocity(tracetemp,params);
        indhere = 1:size(vel_temp,1);
vel_temp(:,mm) = veloc(indhere);
std_temp(:,mm) = stand(indhere);

accel_temp(:,mm) = accelout(indhere);
av_std_accel(:,mm) = accelstd(indhere);

    end
    vel_comps(:,jj) = nanmean(vel_temp,2);
    av_vel(:,jj) =  nanmean(vel_temp,2).^2;
    
      std_comps(:,jj) = nanmean(std_temp,2);
    av_std(:,jj) =  nanmean(std_temp,2).^2;
    
    %% accel
 accel_comps(:,jj) = nanmean(accel_temp,2);
    av_accel(:,jj) =  nanmean(accel_temp,2).^2;
    
      std_comps_accel(:,jj) = nanmean(av_std_accel,2);
    av_std_accel(:,jj) =  nanmean(av_std_accel,2).^2;
end
av_std = sqrt(nanmean(av_std,2));
av_vel = sqrt(nanmean(av_vel,2));

%%accel
av_std_accel = sqrt(nanmean(av_std_accel,2));
av_accel = sqrt(nanmean(av_accel,2));
end
