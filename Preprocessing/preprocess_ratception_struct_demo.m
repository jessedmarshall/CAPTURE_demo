function [ratception_struct,analysisstruct_app] = preprocess_ratception_struct_demo(aa,preprocessing_parameters,params)

%% load in params for back compatability

%links_here =params.links_here;

%mattag = 'JDM37_dannce';


repfactor = params.repfactor ;
%loading


%% preprocess and generate a mocapstruct
ratception_struct = [];

 markernames = fieldnames(aa.predictions);

%% replicate
ratception_struct.fps = 300;

%% don't have rest
ratception_struct.markers_preproc = aa.predictions;

ratception_struct.markernames = markernames;
%test = gtruthmocap.predictions;

ratception_struct.markers = ratception_struct.markers_preproc;
ratception_structtemp = ratception_struct;
chunksize = 10^5;
totalsize = size(ratception_struct.markers_preproc.SpineF,1);
ratception_struct_full =[];
filelength=0;

for rk =1:ceil(totalsize/chunksize)
ratception_struct_temp = ratception_structtemp;
for ll = 1:numel(markernames)
ratception_struct_temp.markers.(markernames{ll}) = ...
    repelem(ratception_struct_temp.markers.(markernames{ll})(1+chunksize*(rk-1):min(chunksize*rk,totalsize),:),repfactor,1);
end
ratception_struct_temp.markers_preproc = ratception_struct_temp.markers;

ratception_struct_temppreproc = compute_preprocessed_mocapstruct(ratception_struct_temp,preprocessing_parameters);


if (rk==1)
    ratception_struct_full = ratception_struct_temppreproc;
    filelength = size(ratception_struct_temppreproc.aligned_mean_position,1);

else
   for ll = 1:numel(markernames)
ratception_struct_full.markers_preproc.(markernames{ll}) = cat(1,ratception_struct_full.markers_preproc.(markernames{ll}),...
     repelem(ratception_struct_temppreproc.markers_preproc.(markernames{ll}),1,1));
 ratception_struct_full.markers_aligned_preproc.(markernames{ll}) = cat(1,ratception_struct_full.markers_aligned_preproc.(markernames{ll}),...
     repelem(ratception_struct_temppreproc.markers_aligned_preproc.(markernames{ll}),1,1));
 ratception_struct_full.bad_frames_agg{ll} = cat(2, ratception_struct_full.bad_frames_agg{ll},...
     filelength+ratception_struct_temppreproc.bad_frames_agg{ll});
end 
    
    ratception_struct_full.move_frames= cat(2,ratception_struct_full.move_frames,filelength+...
        repelem(ratception_struct_temppreproc.move_frames,1,1));% ratception_struct.markers_aligned_preproc.(markernames{ll})+10;
ratception_struct_full.rest_frames= cat(2,ratception_struct_full.rest_frames,filelength+...
    repelem(ratception_struct_temppreproc.rest_frames,1,1));% ratception_struct.markers_aligned_preproc.(markernames{ll})+10;
ratception_struct_full.aligned_rotation_matrix= cat(3,ratception_struct_full.aligned_rotation_matrix,...
    repelem(ratception_struct_temppreproc.aligned_rotation_matrix,1,1,1));% ratception_struct.markers_aligned_preproc.(markernames{ll})+10;
ratception_struct_full.aligned_mean_position= cat(1,ratception_struct_full.aligned_mean_position,repelem(ratception_struct_temppreproc.aligned_mean_position,1,1));% ratception_struct.markers_aligned_preproc.(markernames{ll})+10;
        filelength = filelength+size(ratception_struct_temppreproc.aligned_mean_position,1);

end
end
ratception_struct=ratception_struct_full;


end