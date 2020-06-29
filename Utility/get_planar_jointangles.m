   
      function [angle_inplane,angle_outofplane] =  get_planar_jointangles(mocapstruct,angle_set_plane,angle_set_compare)
          fname_aligned = fieldnames(mocapstruct.markers_aligned_preproc);
          %% get vectors and normalize
          for ii=1:2
          if ~strcmp(angle_set_plane{ii}{1},'zvector')
       vector_plane{ii} =mocapstruct.markers_aligned_preproc.(angle_set_plane{ii}{1})(:,:)-...
           mocapstruct.markers_aligned_preproc.(angle_set_plane{ii}{2})(:,:);
          else
                vector_plane{ii} = ones(size(mocapstruct.markers_aligned_preproc.(fname_aligned{1})(:,:)));
                 vector_plane{ii}(:,[1,2]) = 0;
          end
          end
            vector_examine =mocapstruct.markers_aligned_preproc.(angle_set_compare{1})(:,:)-...
           mocapstruct.markers_aligned_preproc.(angle_set_compare{2})(:,:);
          
          for ii=1:2
               vector_plane{ii} =  vector_plane{ii}./sqrt(nansum( vector_plane{ii}.^2,2));
          end
          %get orthonormal basis for in and out of plane compoenents 
          vector_plane{2} = vector_plane{2}-dot(vector_plane{2}',vector_plane{1}')'.*vector_plane{2};
          vector_plane{2} =  vector_plane{2}./sqrt(sum( vector_plane{2}.^2,2));
     vector_outofplane = cross(vector_plane{1},vector_plane{2});
     vector_outofplane = vector_outofplane./sqrt(nansum(vector_outofplane.^2,2)); %probably not necessary
%get vector to examine
           vector_examine=  vector_examine./sqrt(nansum( vector_examine.^2,2));
    
           %% represent the in plane as a combination of the basis vectors
     vector_examine_inplane =  dot(vector_examine',vector_plane{1}')'.*vector_plane{1}+...
         dot(vector_examine',vector_plane{2}')'.*vector_plane{2};
          angle_inplane = real(acosd(dot(vector_examine_inplane',vector_plane{1}')'./sqrt(sum( vector_examine_inplane.^2,2))));
     
      %% represent the in plane as a combination of the basis vectors
     vector_examine_outofplane =  dot(vector_examine',vector_plane{1}')'.*vector_plane{1}+...
         dot(vector_examine',vector_outofplane')'.*vector_outofplane;
          angle_outofplane = real(acosd(dot(vector_examine_outofplane',vector_plane{1}')'./sqrt(nansum( vector_examine_outofplane.^2,2))));
     figure(68)
plot(angle_inplane(mocapstruct.move_frames_fast))
hold on
plot(angle_outofplane(mocapstruct.move_frames_fast),'r')
hold off
      end