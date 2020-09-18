
function M = plot_multi_clusters_repro_cont_bird(mocapstruct,ind_cell,nframes,v,fighand,titlein)

%num2str(cluster_numbers');
if nargin <5
    fighand= gcf;%   figure(fighand)
else
    fighand =   figure(388);
end
set(fighand,'Color','k')
set(fighand,'Position',[100 100 1100 1100])
nrows = 2;
ncols = 3;
nreps =5;
num_vids = min(numel(ind_cell),nrows*ncols);
maxreps = min(nreps,ceil(numel(ind_cell)./(nrows*ncols)));

for ll = 1:num_vids  %numel(cluster_numbers)
    %  mov_ind = cluster_numbers(ll);
    %     h{ll}=subplot_tight(nrows,ncols, ll);
    %  ind_cell{ll} = cat(1,ind_cell{ll},max(ind_cell{ll})+(5:5:20)');
end
for zz = 1:maxreps
    M_indiv = cell(1,nrows*ncols);
    for ll = 1:num_vids  %numel(cluster_numbers)
        %       subplot_tight(nrows,ncols, ll);
        viduse = (zz-1)*ncols*nrows+ll;
        % if this is a valid frame
        if viduse <=numel(ind_cell)
            
            movie_size = numel(ind_cell{viduse});
            %if any frames at all
            if movie_size>0
                
                if movie_size>=nframes
                    frames_plot_here = ind_cell{viduse}(1:nframes);
                else
                    frames_plot_here = cat(1,ind_cell{viduse},ind_cell{viduse}*ones(1,nframes-movie_size));
                end
                clf;
                %get the individual movies
                M_indiv{ll} =  plot_frame_dannce_reprojection_multi_bird(mocapstruct,frames_plot_here,...
                    mocapstruct.basefolder,mocapstruct.camuse);
                axis off
                %if empty
            else
                set(gca,'Color','k')
            end
            
        end
    end
    
    badinds = find(cellfun(@numel,M_indiv)==0);
    for mm=badinds
        M_indiv{mm} = M_indiv{1};
    end
    
    %% concatenate -- currently hard coded for 6
    M_agg = M_indiv{1}(1);
    for lk=1:nframes%numel(M_agg)
        M_agg(lk).cdata = cat(1,cat(2,M_indiv{1}(lk).cdata,M_indiv{2}(lk).cdata,M_indiv{3}(lk).cdata),...
            cat(2,M_indiv{4}(lk).cdata,M_indiv{5}(lk).cdata,M_indiv{6}(lk).cdata));
    end
    
    doubleit=1;
    if doubleit
            for lk=1:nframes%numel(M_agg)
  M_agg(lk+nframes).cdata  = M_agg(lk).cdata;
            end
    end
    
   M=[];
    %M_agg = cat(1,M_indiv{1},M_indiv{2},M_indiv{3}
    %concatenate
    
    % add title
    if nargin>5
             ntitle(titlein  ,'Color','w','FontSize',26,'location','north' )
    end
    M = M_agg;
    writeVideo(v,M_agg)
   % M_agg =[];
    
    
end
clf;
end

%v = VideoWriter(strcat(savedirectory_subcluster,'aggregate_movie_2'),'MPEG-4');
%     open(v)
%     writeVideo(v, M_here)
%     close(v)