function plot_frame_sub(mocapstruct,frame_inds,h,nmarkers)
%matlab_fr = 10;
%h=figure(370)
% frame_inds = time_ordering_fulltrace{jjj}(1:matlab_fr:min(frames_use,numel(time_ordering_fulltrace{jjj})));

%M{jjj} = movie;
%Mhere = movie;
frame_last = 0;

if nargin<4
    nmarkers = numel(mocapstruct.markernames);
end

marker_plot = ones(1,numel(mocapstruct.markernames));


%% initialize the figure
%    set(h,'Color','k')
%     plot3( squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,1)),...
%         squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,2)),...
%         squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,3)),'o','Color',mocapstruct.markercolor{1},'MarkerFaceColor',mocapstruct.markercolor{1},'MarkerSize',6)
%    
    
xx = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,1));
yy = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,2));
zz = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,3));
line(xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{1},'MarkerFaceColor',mocapstruct.markercolor{1},'MarkerSize',6);

    
    ax = gca;
    axis(ax,'manual')
    set(gca,'Color','k')
  %  grid on;
    set(gca,'Xcolor',[1 1 1 ]);
    set(gca,'Ycolor',[1 1 1]);
    set(gca,'Zcolor',[1 1 1]);
    
%     zlim([-110 120])
%     xlim([-120 120])
%     ylim([-120 120])
%         


%     xlim([-135 135])
%     ylim([-135 135])
% % 
    zlim([-110 110])
    xlim([-135 135])
    ylim([-135 135])
    set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])
        %view([-22, 12]);
        view([162, 7]);

        base_time = datenum(mocapstruct.mocapfiletimes{1});

        
  
       % datestr(datenum(mocapstruct_social.mocapfiletimes{1})+300./(300*60*60*24))
for lk = reshape(frame_inds,1,[])%1:10:10000
     cla;
    
  %  set(handles.t4,'String',num2str(lk));
    
msize = 6;
    
    ind_to_plot = lk;

    %% Plot markers that are tracked in the frame
  set(gca,'Nextplot','ReplaceChildren');
    handles_here = cell(1,numel(mocapstruct.markernames));
    for jj = 1:nmarkers
        % don't plot markers that drop out
        if ~isnan(sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2))
        if (~sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2) == 0)
%             handles_here{jj} = plot3( squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,1)),...
%                 squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,2)),...
%                 squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,3)),...
%'o','Color',mocapstruct.markercolor{jj},'MarkerFaceColor',mocapstruct.markercolor{jj},'MarkerSize',8);
%         
%             xx = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,1));
%                 yy = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,2));
%                 zz = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,3));
%                 handles_here{jj} = line(xx,yy,zz,'Linestyle','-','Color',mocapstruct.markercolor{jj},'Marker','o','LineWidth',2,'Color',mocapstruct.markercolor{jj},'MarkerFaceColor',mocapstruct.markercolor{jj},'MarkerSize',6);
%                     line(xx,yy,zz,'Linestyle','-','Color',mocapstruct.markercolor{jj},'LineWidth',2);

                           hold on

            
           
            marker_plot(jj) = 1;
        else
            marker_plot(jj) = 0;
        end
        end
    end
    
    msize = 6;
    %% plot the links between markers
    for mm = numel(mocapstruct.links):-1:1
        if numel(   mocapstruct.links{mm})
        if   mocapstruct.links{mm}(1)<=nmarkers && mocapstruct.links{mm}(2)<=nmarkers
        if (ismember(mocapstruct.links{mm}(1),1:numel(mocapstruct.markernames)) && ismember(mocapstruct.links{mm}(2),1:numel(mocapstruct.markernames)))
            if (marker_plot(mocapstruct.links{mm}(1)) == 1 && marker_plot(mocapstruct.links{mm}(2)) == 1)
                
                
                  xx = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,1)) ...
                    squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,1)) ];
                yy = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,2)) ...
                    squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,2))];
                zz = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,3)) ...
                    squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,3))];
               % line(xx,yy,zz,'Color',mocapstruct.markercolor{mocapstruct.links{mm}(1)},'LineWidth',2);
               %  line(xx,yy,zz,'Color',mocapstruct.markercolor{mocapstruct.links{mm}(1)},'LineWidth',2,'Marker','o','LineWidth',2,'Color',...
                %    mocapstruct.markercolor{mocapstruct.links{mm}(1)},'MarkerFaceColor',mocapstruct.markercolor{mocapstruct.links{mm}(1)},'MarkerSize',6)
                      line(xx,yy,zz,'Color',mocapstruct.markercolor{mocapstruct.links{mm}(1)},'LineWidth',2,'LineWidth',2,'Color',...
                    mocapstruct.markercolor{mocapstruct.links{mm}(1)})%,'MarkerFaceColor',mocapstruct.markercolor{mocapstruct.links{mm}(1)},'MarkerSize',6)
%                 
                
%                 plot3( [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,1)) ...
%                     squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,1)) ],...
%                     [ squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,2))...
%                     squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,2))],...
%                     [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,3))...
%                     squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,3))],...
%                     'Color',mocapstruct.markercolor{mocapstruct.links{mm}(1)},'Linewidth',3);
             end
        end
        end
        end
    end

        handles_here = cell(1,numel(mocapstruct.markernames));

    for jj = 1:nmarkers
        % don't plot markers that drop out
        if ~isnan(sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2))
        if (~sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2) == 0)
             handles_here{jj} =     plot3( squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,1)),...
                squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,2)),...
                squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,3)),...
'o','Color',mocapstruct.markercolor{jj},'MarkerFaceColor',mocapstruct.markercolor{jj},'MarkerSize',msize);
        end
        end
    end


    
        %title(strcat('  Frame: ' ,datestr(base_time+lk./(mocapstruct.fps*60*60*24))),'Color','w')
    
    
   %drawnow
    %hold off
  
    frame_last = lk;
    set(gca,'YTickLabels','')
    axis off
    
      %  M(find(frame_inds == lk)) =  getframe(gcf);

    %clf
end
%set(gca,'Nextplot','add');

end