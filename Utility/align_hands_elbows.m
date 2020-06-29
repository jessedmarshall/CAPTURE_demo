function [markers, markers_aligned,cluster_mean_array,global_rotmatrix] = align_hands_elbows(markers,fps)


markers_aligned = struct();
markernames = fieldnames(markers);

%% remove high frequency noise (>10 Hz) in the spinal trace
% 
% dH = designfilt('lowpassiir', 'FilterOrder', 3, 'HalfPowerFrequency', 10/(fps/2), ...
%     'DesignMethod', 'butter');
% [f1,f2] = tf(dH);
% 
% marker_F_lowpass =zeros(size(markers.SpineF));
% marker_M_lowpass =zeros(size(markers.SpineM));
% 
% %% shouldn't be any nans, but if there are they screw up filtfilt
% markers.SpineF(isnan(markers.SpineF)) = 0;
% markers.SpineM(isnan(markers.SpineM)) = 0;
% 
% 
% for mk = 1:3
% marker_F_lowpass(:,mk) = filtfilt(f1,f2,...
%             markers.SpineF(:,mk));
%         
% marker_M_lowpass(:,mk) = filtfilt(f1,f2,...
%             markers.SpineM(:,mk));
% end

marker_F_lowpass = markers.SpineF;
marker_M_lowpass = markers.SpineM;

%% get the marker rotation matrix
fprintf('rotating markers')
rotangle = atan2(-(marker_F_lowpass(:,2)-marker_M_lowpass(:,2)),...
    (marker_F_lowpass(:,1)-marker_M_lowpass(:,1)));
global_rotmatrix = zeros(2,2,numel(rotangle));
global_rotmatrix(1,1,:) = cos(rotangle);
global_rotmatrix(1,2,:) = -sin(rotangle);
global_rotmatrix(2,1,:) = sin(rotangle);
global_rotmatrix(2,2,:) = cos(rotangle);
cluster_mean_array =  (markers.SpineM(: ,:));


%% subtract the mean and rotate
for mm = 1:numel(markernames)
    markers_aligned.(markernames{mm}) = markers.(markernames{mm});
end


for mm = (1:numel(markernames))
    markers_aligned.(markernames{mm}) =  bsxfun(@minus,...
        markers_aligned.(markernames{mm}), cluster_mean_array(:,:));
end


for mm = (1:numel(markernames))
    markers_aligned.(markernames{mm})(:,1:2) =  ...
        squeeze(mtimesx(global_rotmatrix(:,:,:),(reshape(markers_aligned.(markernames{mm})(:,1:2)',2,1,[]))...
        ))';
end


for mm = (1:numel(markernames))
    markers_aligned.(markernames{mm})(find(isnan(markers_aligned.(markernames{mm})))) =  0;
        end
%look at at most 50 clusters
%
% figure(44)
% plot(markers_aligned.ArmL(:,3),'r')
% hold on
% plot(markers_aligned.ElbowL(:,3),'b')
% hold off
do_swap = 0;
if do_swap
if isfield(markers_aligned,'ArmL')
vec1 = marker_F_lowpass(:,[1,3])-marker_M_lowpass(:,[1,3]);
vec2 = markers_aligned.ArmL(:,[1,3])-markers_aligned.ElbowL(:,[1,3]);
norm1 = sqrt(sum(vec1.^2,2));
norm2 = sqrt(sum(vec2.^2,2));

dp = acosd(dot(vec1',vec2')./(norm1.*norm2)');


frameswap = find(dp>110);
temp_elbow = markers_aligned.ElbowL;
markers_aligned.ElbowL(frameswap,:) = markers_aligned.ArmL(frameswap,:);
markers_aligned.ArmL(frameswap,:) = temp_elbow(frameswap,:);
clear temp_elbow

vec1 = marker_F_lowpass(:,[1,3])-marker_M_lowpass(:,[1,3]);
vec2 = markers_aligned.ArmL(:,[1,3])-markers_aligned.ElbowL(:,[1,3]);
norm1 = sqrt(sum(vec1.^2,2));
norm2 = sqrt(sum(vec2.^2,2));

dp_new = acosd(dot(vec1',vec2')./(norm1.*norm2)');
%
% bad_frames_elbowspine = unique(cat(2,bad_frames_agg{4},bad_frames_agg{5},bad_frames_agg{11},bad_frames_agg{12}));
% goodframeshere =setxor(1:size(vec1,1),bad_frames_elbowspine);

%                                 figure(33)
%
%                                 plot(dp(goodframeshere),'r')
%                                 hold on
%                                 plot(dp_new(goodframeshere),'b')
%                                 hold off


markers_arml_pre = markers.ArmL;
markers_elbow_pre=markers.ElbowL;
markers.ArmL(frameswap,:) = markers.ElbowL(frameswap,:);
markers.ElbowL(frameswap,:) = markers_arml_pre(frameswap,:);
end
end
end