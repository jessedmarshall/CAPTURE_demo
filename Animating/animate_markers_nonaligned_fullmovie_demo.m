function M = animate_markers_nonaligned_fullmovie_demo(mocapstruct,frame_inds,fighand,axisparams)
%matlab_fr = 10;
if nargin<3 
elseif isempty(fighand)
    h=figure(370)
else
    h=fighand;
end

%% the default case
 [maxval] = max(mocapstruct.markers_preproc.SpineM,[],1);
  [minval] =   min(mocapstruct.markers_preproc.SpineM,[],1);
buffactor_axis = 1.3; %1.1
buffactor_arena = 1.02;


    th = 0:pi/100:2*pi;
    xcent = (maxval(1)+minval(1))./2;
      ycent = (maxval(1)+minval(1))./2;

  

      
if nargin<4 || isempty(axisparams)
          zlimvals = [-50 250];
              xlimvals = [-(304*buffactor_axis- xcent) (304*buffactor_axis+ xcent)];
      ylimvals = [-(304*buffactor_axis- ycent) (304*buffactor_axis+ ycent)];
else
    zlimvals = axisparams.zlim ;
        xlimvals = axisparams.xlim ;
    ylimvals = axisparams.ylim ;

end


frame_last = 0;

marker_plot = ones(1,numel(mocapstruct.markernames));


%% initialize the figure
    set(h,'Color','k')

xx = squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{1})(1,1));
yy = squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{1})(1,2));
zz = squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{1})(1,3));
handle_base = line(xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{1},'MarkerFaceColor',mocapstruct.markercolor{1},'MarkerSize',6);

    
    ax = gca;
    axis(ax,'manual')
    set(gca,'Color','k')
    grid off;
    axis off
    set(gca,'Xcolor',[1 1 1 ]);
    set(gca,'Ycolor',[1 1 1]);
    set(gca,'Zcolor',[1 1 1]);
    
  %  zlim([200 350])
   % xlim([-340 340])
   % ylim([-340 340])
    set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])
        view([-22, 12]);
        
        
        
  
xunit = 304*buffactor_arena * cos(th) + xcent;
yunit = 304*buffactor_arena * sin(th) + ycent;
zunit = zeros(1,numel(xunit));
plot3(xunit,yunit,zunit,'w','linewidth',3);

 zlim(zlimvals)
    xlim(xlimvals)
    ylim(ylimvals)
    set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])

axis tight
        


    set(gca,'Xcolor',[1 1 1 ]);
    set(gca,'Ycolor',[1 1 1]);
    set(gca,'Zcolor',[1 1 1]);
    

   ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width 1.1*ax_height];
        view([162, 7]);

for lk = reshape(frame_inds,1,[])%1:10:10000
    cla;
   fprintf('frame: %f \n', lk)

    plot3(xunit,yunit,zunit,'w','linewidth',3);

    ind_to_plot = lk;

    %% Plot markers that are tracked in the frame
  set(gca,'Nextplot','ReplaceChildren');
    handles_here = cell(1,numel(mocapstruct.markernames));
    for jj = 1:numel(mocapstruct.markernames)
        % don't plot markers that drop out
        if ~isnan(sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2))
        if (~sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2) == 0)      
            xx = squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,1));
                yy = squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,2));
                zz = squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,3));
                handles_here{jj} = line(xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{jj},'MarkerFaceColor',mocapstruct.markercolor{jj},'MarkerSize',6);
                
                
            
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
                  
                  xx = [squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,1)) ...
                      squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,1)) ];
                  yy = [squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,2)) ...
                      squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,2))];
                  zz = [squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,3)) ...
                      squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,3))];
                  line(xx,yy,zz,'Color',mocapstruct.markercolor{mocapstruct.links{mm}(1)},'LineWidth',2);
              end
              
          end
      end
  end
  
 zlim(zlimvals)
    xlim(xlimvals)
    ylim(ylimvals)
    set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])

   
axis off

    
    drawnow 
    hold off
  
    frame_last = lk;
    
        M(find(frame_inds == lk)) =  getframe(gcf);

    %clf
end

end