function bad_frames_agg = getbadframes(marker_velocity,marker_position,fps,preprocessing_parameters)
num_markers = size(marker_velocity,1);


  bad_frame_velthresh = preprocessing_parameters.bad_frame_vel_thresh;
      bad_frame_surround_flag =   preprocessing_parameters.bad_frame_surround_flag;
  bad_frame_surround=  preprocessing_parameters.bad_frame_surround_number;

  
for mm = 1:num_markers
    bad_1 = find(squeeze(marker_velocity(mm, :,4))>bad_frame_velthresh);
    [i3] = ind2sub(size(squeeze(marker_velocity(mm, :,4))),bad_1);
    clear bad_1
    
    % find zeros - no marker found
    bad_2 = find(squeeze(sum(marker_position(mm,:, :),3))==0);
    [i4] = ind2sub(size(squeeze(sum(marker_position(mm,:, :),3))),bad_2);
    clear bad_2
    
    % nans can also be no marker found after some processing
       bad_3 = find(squeeze(sum(isnan(marker_position(mm,:, :)),3)));
    [i5] = ind2sub(size(squeeze(sum(~isnan(marker_position(mm,:, :)),3))),bad_3);
    clear bad_3
    
    if numel(i4)
        bad_frames  = (cat(2,(i3),(i4)));
    else
        bad_frames = i3;
    end
    clear i3 i4
    
     if numel(i5)
        bad_frames = cat(2,bad_frames,i5);
     end
    bad_frames = unique((bad_frames));
    numframes = size(marker_velocity,2);
    gap_base = zeros(1,numframes );
    gap_base(bad_frames) = 1;
    
    if bad_frame_surround_flag
    gap_base_conv = conv(gap_base,ones(1,bad_frame_surround)./bad_frame_surround,'same');
    %gap_base_conv(gap_base_conv>0) = 1;
    bad_frames = find(gap_base_conv >0);
    else
      bad_frames = find(gap_base>0);   
    end
    bad_frames = reshape(bad_frames,1,[]);
    bad_frames = unique(bad_frames);
    bad_frames(bad_frames<1)=1;
    bad_frames(bad_frames>numframes ) = numframes ;
    
    bad_frames_agg{mm} = bad_frames;
end

end