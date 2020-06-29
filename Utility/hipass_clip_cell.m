function [marker_clipped_return,clipped_index] = hipass_clip_cell(markers,fragment_frames,params,subsetnames)

dHipass = designfilt('highpassiir', 'FilterOrder', 3, 'HalfPowerFrequency', 0.3/(params.fps/2), ...
    'DesignMethod', 'butter');
[f1_hipass,f2_hipass] = tf(dHipass);

marker_clipped = markers;

inst_label = zeros(1,size(markers,1));
inst_label(fragment_frames) = 1;
clipped_index = 1:numel(inst_label);
pixellist = bwconncomp(inst_label);

%% interpolation utility
  interp_length = 5;
            interp_x = [1:interp_length (2*interp_length+1):(3*interp_length)];
            interp_desired = (interp_length+1):(2*interp_length);

for mk = 1:size(markers,2)
    fprintf('Clipping marker %s \n',num2str(mk))
    for lk =1:pixellist.NumObjects
        %params.min_frames;

        if (numel(pixellist.PixelIdxList{lk})>=1) % at least 200 ms consecutive good tracking to be used -- removed to facilitate other comparisons
                 if (numel(pixellist.PixelIdxList{lk})>=10)
        frame_fragment = filtfilt(f1_hipass,f2_hipass,marker_clipped((pixellist.PixelIdxList{lk}),mk)-...
                mean(marker_clipped((pixellist.PixelIdxList{lk}),mk)));
          
            
            interp_y = cat(2,zeros(1,interp_length),frame_fragment(interp_desired)');
            interp_output = spline(interp_x,interp_y,interp_desired);
            frame_fragment(1:interp_length) = interp_output;
            
            fragment_length = length(frame_fragment);
             interp_y_end= cat(2,frame_fragment((fragment_length-2*interp_length+1):(fragment_length-interp_length))',zeros(1,interp_length));
                         interp_output_end = spline(interp_x,interp_y_end,interp_desired);
 frame_fragment((fragment_length-interp_length+1):(fragment_length)) = interp_output_end;
            
            
                        marker_clipped((pixellist.PixelIdxList{lk}),mk) = frame_fragment;
                 else 
                     %don't hipass if not enough frames
                      marker_clipped((pixellist.PixelIdxList{lk}),mk) = marker_clipped((pixellist.PixelIdxList{lk}),mk)-...
                mean(marker_clipped((pixellist.PixelIdxList{lk}),mk));
            
            
                 end
            
        
        else
            inst_label(pixellist.PixelIdxList{lk}) = 0;
        end
        
    end
    
    
    %% also remove the bad frames
end

%% interpolate across the bad frames

marker_clipped_return = marker_clipped(find(inst_label==1),:);
%for mk =  1:size(markers,2)
      marker_clipped(find(inst_label == 0),:) = [];
%end
    clipped_index(inst_label==0) = [];


end