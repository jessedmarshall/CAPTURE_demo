function M = animate_markers_aligned_fullmovie_demo(mocapstruct,frame_inds,fighand)
%matlab_fr = 10;
if nargin<3
h=figure(370)
else
    h=fighand;
end

frame_last = 0;

marker_plot = ones(1,numel(mocapstruct.markernames));


%% initialize the figure
    set(h,'Color','k')
    
xx = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,1));
yy = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,2));
zz = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,3));
handle_base = line(xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{1},'MarkerFaceColor',mocapstruct.markercolor{1},'MarkerSize',6);

    
    ax = gca;
    axis(ax,'manual')
    set(gca,'Color','k')
    grid on;
    set(gca,'Xcolor',[1 1 1 ]);
    set(gca,'Ycolor',[1 1 1]);
    set(gca,'Zcolor',[1 1 1]);
    
    zlim([-110 170])
    xlim([-140 140])
    ylim([-140 140])
    
%     zlim([-210 270])
%     xlim([-240 240])
%     ylim([-240 240])
%     
 
    set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])
        view([-22, 12]);


for lk = reshape(frame_inds,1,[])%1:10:10000
    fprintf('frame %f \n',lk);
    cla;
    
    
 %mocapstruct.links{20} = [];
 %       mocapstruct.links{22} = [];

    ind_to_plot = lk;

    %% Plot markers that are tracked in the frame
  set(gca,'Nextplot','ReplaceChildren');
    handles_here = cell(1,numel(mocapstruct.markernames));
    for jj = 1:numel(mocapstruct.markernames)
        % don't plot markers that drop out
        if ~isnan(sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2))
        if (~sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2) == 0)      
            xx = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,1));
                yy = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,2));
                zz = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,3));
                handles_here{jj} = line(xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{jj},'MarkerFaceColor',mocapstruct.markercolor{jj},'MarkerSize',9);
                
                
            
            hold on
            marker_plot(jj) = 1;
        else
            marker_plot(jj) = 0;
        end
        
        end
    end
   
  %% plot the links between markers
  for mm = 1:numel(mocapstruct.links)
      if numel(mocapstruct.links{mm})
          if (ismember(mocapstruct.links{mm}(1),1:numel(mocapstruct.markernames)) && ismember(mocapstruct.links{mm}(2),1:numel(mocapstruct.markernames)))
              if (marker_plot(mocapstruct.links{mm}(1)) == 1 && marker_plot(mocapstruct.links{mm}(2)) == 1)
                  
                  xx = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,1)) ...
                      squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,1)) ];
                  yy = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,2)) ...
                      squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,2))];
                  zz = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,3)) ...
                      squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,3))];
                  line(xx,yy,zz,'Color',mocapstruct.markercolor{mocapstruct.links{mm}(1)},'LineWidth',3);
              end
              
          end
      end
  end

    drawnow 
    hold off
  
    frame_last = lk;
    
        M(find(frame_inds == lk)) =  getframe(gcf);

    %clf
end
%set(gca,'Nextplot','add');

end