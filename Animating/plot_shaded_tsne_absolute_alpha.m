function h2 = plot_shaded_tsne_absolute_alpha(analysisstruct,shading,names,iter,figinput,sortedind)
 % subplot(1,2,1)
% for nn = 1:numel(conditions)
%       plot(0,0,'o','MarkerSize',8,'Color',colors(nn,:),'MarkerFaceColor',colors(nn,:))
%     hold on
% end
if nargin<5
h2=figure(606+iter)
else
    h2 = figinput;
end
if nargin<6
    sortedind=1
end


shading_values = 100.*(shading-[min(shading)+0.001])./(max(shading)-[min(shading)+0.001]);
shading_values(shading_values<0) = 0;
if find(isnan(shading_values))
    fprintf('shading is nan? \n')
end
%shading(isnan(shading)) = 0.25;
colors = cat(1,0.95*ones(3,3),othercolor('PuRd9',101));
% numcolors=101-(10+1);
% redmap = cbrewer('seq','Reds',numcolors);
% bluemap = flipud(cbrewer('seq','Blues',numcolors));
% halfcolor = floor(numcolors./2);
% colormaphere = cat(1,bluemap(1:halfcolor,:),ones(10,3),redmap(halfcolor+1:end,:));
% colors = cat(1,0.95*ones(1,3),colormaphere);


set(h2,'Color','w')
     if sortedind
      nnn = analysisstruct.sorted_watershed;
nnn(nnn>0) = 1;
B = bwboundaries((flipud(nnn)));
hold on
for kk = 1:numel(B)
        if kk<=numel(shading_values)
    kkhere = (analysisstruct.sorted_clust_ind(kk));
                %if numel(find(analysisstruct.annot_reordered{end,end}==find(analysisstruct.sorted_clust_ind==kk)))>plotthresh
     fill(analysisstruct.xx(B{kkhere}(:,2)),analysisstruct.yy(numel(analysisstruct.yy)-B{kkhere}(:,1)),colors(1+floor(shading_values(kk)),:),...
        'EdgeColor',colors(1+floor(shading_values(kk)),:)); %,'none'
               % end
    end
end
hold off
     else
     nnn = analysisstruct.unsorted_watershed;
     fprintf('running unsorted \n')
nnn(nnn>0) = 1;
B = bwboundaries((flipud(nnn)));
hold on
for kk = 1:numel(B)
        if kk<=numel(shading_values)
                %if numel(find(analysisstruct.annot_reordered{end,end}==find(analysisstruct.sorted_clust_ind==kk)))>plotthresh
    fill(analysisstruct.xx(B{kk}(:,2)),analysisstruct.yy(numel(analysisstruct.yy)-B{kk}(:,1)),colors(1+floor(shading_values(kk)),:),...
        'EdgeColor',colors(1+floor(shading_values(kk)),:),'Linewidth',2); %,'none'
    %hhere.FaceColor = colors(1+floor(shading_values(kk)),:);
               % end
        end
end
hold off     
     end
     axis equal
% 
% for ll = (1:min(analysisstruct.density_objects,numel(shading_values)))    
%  %   clust_reordered = find(analysisstruct.sorted_clust_ind==ll);
%   %  if numel(clust_reordered)
%   if sortedind
%     plot(analysisstruct.zValues(find((analysisstruct.annot_reordered{end,end})==ll),1),...
%        analysisstruct.zValues(find(analysisstruct.annot_reordered{end,end}==ll),2),'o','MarkerSize',2,'MarkerFaceColor',colors(1+floor(shading_values(ll)),:),'MarkerEdgeColor','none');
%   else
%      plot(analysisstruct.zValues(find((analysisstruct.annotation_vec{end,end})==ll),1),...
%        analysisstruct.zValues(find(analysisstruct.annotation_vec{end,end}==ll),2),'o','MarkerSize',2,'MarkerFaceColor',colors(1+floor(shading_values(ll)),:),'MarkerEdgeColor','none');  
%   end
%    
%    % end 
%    hold on
% end
hold off
xlabel('Tsne 1')
ylabel('Tsne 2')
set(gca,'fontsize',18)
box off
colormap(colors)
colorbar
cbh=colorbar;
set(cbh,'YTick',[0 0.5 1.0],'TickLabels',{num2str(round(min(shading),1)) num2str(round(min(shading)+(max(shading)-min(shading))./2,1)) num2str(round(max(shading),1))})

if numel(names)
%% plot names as well
hold on
s2 = regionprops(analysisstruct.unsorted_watershed, 'Centroid');
for k = 1:numel(analysisstruct.sorted_clust_ind')
    % get the ind of the sorted cluster
    c = s2((k)).Centroid;
    text(analysisstruct.xx(floor(c(1))), analysisstruct.yy(floor(c(2))), num2str((names(k))), ... %sprintf('%d', integer(ind)),
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','Color','Red','FontWeight','Bold');
end
hold off
end

if isfield(analysisstruct,'coarse_borders')
hold on
B = analysisstruct.coarse_borders
for kk = 1:numel(B)
        plot(analysisstruct.xx((B{kk}{1}(:,2))),...
            analysisstruct.yy((B{kk}{1}(:,1))),...
            'Color',[0 0 0],'linewidth',0.01)

end
hold off
end


%if numel(legendinput)
%legend(legendinput)
%end
end