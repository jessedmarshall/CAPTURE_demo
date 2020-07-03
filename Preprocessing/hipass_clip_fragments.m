function [marker_clipped,clipped_index] = hipass_clip_fragments(markers,fragment_frames,params,subsetnames)

dHipass = designfilt('highpassiir', 'FilterOrder', 3, 'HalfPowerFrequency', 0.3/(params.fps/2), ...
    'DesignMethod', 'butter');
[f1_hipass,f2_hipass] = tf(dHipass);

markernames = fieldnames(markers);


if nargin>3
    marker_clipped = struct();
    looplist = subsetnames;
    for jj = subsetnames
    marker_clipped.(markernames{(jj)}) = markers.(markernames{(jj)});
    end
else

marker_clipped = struct();
marker_clipped = markers;
looplist = 1:numel(markernames);
end
%cluster_size_dist = cell(1,number_of_clust(mmm));

inst_label = zeros(1,numel(markers.(markernames{1})(:,1)));
inst_label(fragment_frames) = 1;
clipped_index = 1:numel(inst_label);
pixellist = bwconncomp(inst_label);

%% interpolation utility
  interp_length = 5;
            interp_x = [1:interp_length (2*interp_length+1):(3*interp_length)];
            interp_desired = (interp_length+1):(2*interp_length);

for mk = looplist
    fprintf('Clipping marker %s \n',markernames{mk})
    for lk =1:pixellist.NumObjects
        %params.min_frames;

        if (numel(pixellist.PixelIdxList{lk})>=1) % at least 200 ms consecutive good tracking to be used -- removed to facilitate other comparisons
        for sk = 1:size(marker_clipped.(markernames{mk}),2)
                 if (numel(pixellist.PixelIdxList{lk})>=10)
        frame_fragment = filtfilt(f1_hipass,f2_hipass,marker_clipped.(markernames{mk})(pixellist.PixelIdxList{lk},sk)-...
                mean(marker_clipped.(markernames{mk})(pixellist.PixelIdxList{lk},sk)));
          
            
            interp_y = cat(2,zeros(1,interp_length),frame_fragment(interp_desired)');
            interp_output = spline(interp_x,interp_y,interp_desired);
            frame_fragment(1:interp_length) = interp_output;
            
            fragment_length = length(frame_fragment);
             interp_y_end= cat(2,frame_fragment((fragment_length-2*interp_length+1):(fragment_length-interp_length))',zeros(1,interp_length));
                         interp_output_end = spline(interp_x,interp_y_end,interp_desired);
 frame_fragment((fragment_length-interp_length+1):(fragment_length)) = interp_output_end;
            
            
                        marker_clipped.(markernames{mk})(pixellist.PixelIdxList{lk},sk) = frame_fragment;
                 else 
                     %don't hipass if not enough frames
                      marker_clipped.(markernames{mk})(pixellist.PixelIdxList{lk},sk) = marker_clipped.(markernames{mk})(pixellist.PixelIdxList{lk},sk)-...
                mean(marker_clipped.(markernames{mk})(pixellist.PixelIdxList{lk},sk));
            
            
                 end
            
        end
        else
            inst_label(pixellist.PixelIdxList{lk}) = 0;
        end
        
    end
    
    
    %% also remove the bad frames
end

%% interpolate across the bad frames


for mk = looplist
      marker_clipped.(markernames{mk})(inst_label == 0,:) = [];
end
    clipped_index(inst_label==0) = [];


end