function vectorout_vel= interpolate_ends(vectorout_vel,interp_length)

 % interp_length = params.difforder_movav;
            interp_x = [1:interp_length (2*interp_length+1):(3*interp_length)];
            interp_desired = (interp_length+1):(2*interp_length);
            
            interp_desired(interp_desired>numel(vectorout_vel)) = numel(vectorout_vel);
            
            interp_y = cat(2,zeros(1,interp_length),vectorout_vel(interp_desired)');
            interp_output = spline(interp_x,interp_y,interp_desired);
            vectorout_vel(1:interp_length) = interp_output;
            
            fragment_length = length(vectorout_vel);
            newind_interp = (fragment_length-2*interp_length+1):(fragment_length-interp_length);
            newind_interp(newind_interp>numel(vectorout_vel)) = numel(vectorout_vel);
            newind_interp(newind_interp<1) = 1;
             interp_y_end= cat(2,vectorout_vel(newind_interp)',zeros(1,interp_length));
                         interp_output_end = spline(interp_x,interp_y_end,interp_desired);
                       
                         lastind_interp = (fragment_length-interp_length+1):(fragment_length);
                         lastind_interp(lastind_interp>numel(vectorout_vel)) = numel(vectorout_vel);
                         lastind_interp(lastind_interp<1) =1;
 vectorout_vel(lastind_interp) = interp_output_end;