function [vectorout,vectorstdout,vectorout_accel,vectorstdout_accel]= get_vector_velocity(vectorin,params)

if params.medfiltorder>1
tracesmooth = medfilt1(vectorin,params.medfiltorder);
else
tracesmooth = vectorin;
end
%gfilter = fspecial('gaussian',[50 1], params.gaussorder);
%tracesmoothed = convn(tracesmooth,gfilter,'same');
tracesmoothed = tracesmooth;
endptcond = 'shrink';
vectorout_vel = cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);


%% try interpolating  

            
    vectorout_vel =  interpolate_ends( vectorout_vel,  3*params.difforder_movav);
%vectorout_vel = cat(1,zeros(params.difforder_movav*1,1),vectorout_vel,zeros(params.difforder_movav*1,1));
vectorout = movmean(vectorout_vel,[floor(params.difforder_movav./2) floor(params.difforder_movav./2)] ,'omitnan','Endpoints',endptcond);%cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);
% vectorout = vectorout((params.difforder_movav+1):(end-params.difforder_movav));
% 
% padsize = (numel(vectorin)-numel(vectorout));
% vectorout = cat(1,vectorout,zeros(padsize,1));

    vectorout_vel_big =  interpolate_ends( vectorout_vel,  3*params.difforder_movav);
vectorstdout = movstd(vectorout_vel_big,[floor(params.difforder_movav./2) floor(params.difforder_movav./2)] ,'omitnan','Endpoints',endptcond);%cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);
% vectorstdout = vectorstdout((params.difforder_movav+1):(end-params.difforder_movav));

  %vectorstdout = cat(1,vectorstdout,zeros(padsize,1));

a_diff_order = floor(params.difforder);
vectorout_vel_a = cat(1,zeros((a_diff_order)-2,1),tracesmoothed(a_diff_order:end,1) - tracesmoothed(1:end-a_diff_order+1),1);

vectorout_accel = cat(1,zeros((a_diff_order)-2,1),vectorout_vel_a(a_diff_order:end,1) - vectorout_vel_a(1:end-a_diff_order+1),1);
%vectorout_accel = cat(1,zeros(params.difforder_movav*2,1),vectorout_accel,zeros(params.difforder_movav*2,1));
    vectorout_accel =  interpolate_ends( vectorout_accel,  3*params.difforder_movav);

vectorout_accel = movmean(vectorout_accel,[floor(params.difforder_movav./2) floor(params.difforder_movav./2)] ,'omitnan','Endpoints',endptcond);%cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);
 %vectorout_accel = vectorout_accel((2*params.difforder_movav+1):(end-2*params.difforder_movav));

vectorstdout_accel = movstd(vectorout_accel,[floor(params.difforder_movav./2) floor(params.difforder_movav./2)] ,'omitnan','Endpoints',endptcond);%cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);
 %vectorstdout_accel = vectorstdout_accel((2*params.difforder_movav+1):(end-2*params.difforder_movav));

%padsize = (numel(vectorin)-numel(vectorout_accel));

%vectorout_accel = cat(1,vectorout_accel,zeros(padsize,1));

%padsize = (numel(vectorin)-numel(vectorstdout_accel));
%vectorstdout_accel = cat(1,vectorstdout_accel,zeros(padsize,1));

end