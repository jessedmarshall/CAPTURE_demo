function [density_maps,xx,yy,density_max,RGB] =plotdensitymaps(zValues,num,fighand_in,density_width,density_max,density_res,name_in)

%figure(fighand)

%% density, watershed, images
%density_max = max(zValues{:})*1.3;%;
%density_res = 500;

if nargin <5
for ll =1:numel(zValues)
    density_max_arr(ll,:) = max(zValues{ll});%;
end
density_max = max(max(density_max_arr))*1.3;
end
names_default={'Oranges9','Greens9','Reds9','Greens9','Reds9','Blues9','Greys9','Purples9','Purples9','Greys9'}

if nargin<7
names={'Oranges9','Greens9','Reds9','Greens9','Reds9','Blues9','Greys9','Purples9','Purples9','Greys9'}
if numel(zValues)==1
    names={'Oranges5'};
end

else
    names = {name_in,names_default{:}};
    if numel(zValues)==1
    names={name_in};
end
end
%names={'Oranges9','Greens9','Reds9','Blues9','Purples9','YlGn9','Greys9','PuRd9','BuGn9','PiYG9','PuOr9'}
%     
% 'Blues7','Greys7',...
%     'Greens7','Oranges7','Reds7','Blues7','Purples7','Reds7','Blues7','Greys7'};%PuRd


for ll =1:numel(zValues)

    tsnehere =zValues{ll};

%[xx,yy,density_maps{ll}] = findPointDensity_unnormJDM(tsnehere(:,:),...
 %   density_width,[density_res density_res],[-density_max density_max]);
fprintf('using unnorm')
%findPointDensity_unnormJDM
%%normalized version == looks weird atm
[xx,yy,density_maps{ll}] = findPointDensity(tsnehere(:,:),...
    density_width,[density_res density_res],[-density_max density_max]);
max_density_maps_arr(ll) = max(max(density_maps{ll}));
end
max_density_maps = max(max_density_maps_arr);


for ll =1:numel(zValues)
colorshere{ll} = cat(1,ones(3,3),othercolor(names{ll},256));%parula(256);

if strcmp(class(fighand_in),'matlab.ui.Figure')
    figure(fighand_in)
else
set(gcf,'currentaxes',fighand_in)
end
set(gcf,'Color','w')

 %subplot(1,numel(zValues)+1,ll)
%P = impixel;
aa = imagesc(flipud(density_maps{ll}));
%  cmap = colormap(jet)
cmap=colormap(colorshere{ll});
if numel(density_maps{ll})
 maphere=   density_maps{ll};
 maphere(isnan(maphere)) = 0;
  maphere(isinf(maphere)) = 0;
  
maphere(maphere<10^(-5)) = 0;
 mapherethresh = maphere;
 mapherethresh(maphere<prctile(maphere,1)) = 0;
 maxval = prctile(mapherethresh(mapherethresh>0),99);
 
caxis([0 maxval])
%caxis([min(density_maps{ll}(:)) max_density_maps])
end
c=colorbar
aa.CData = get(aa,'CData');
%cmap=colormap(hh);
c.Label.String = 'Probability Density';
c.Label.FontSize = 18;
c.FontSize = 16;
%title(strrep(conditionnames{ll},'_',''))
box off
cmin = min(aa.CData(:));
cmax = maxval;%prctile(aa.CData(:),99);
m = length(cmap);
index = fix((aa.CData-cmin)/(cmax-cmin)*m)+1; %A
% Then to RGB
RGB{ll} = ind2rgb(index,cmap);

axis off
set(fighand_in,'Color','w')
end

colorshere{ll} = cat(1,ones(3,3),othercolor('Oranges5',256));%parula(256);
plot_axes = 1;
if plot_axes
    colorshere_gray = cat(1,ones(3,3),othercolor('Greys5',256));%parula(256);

figure(479)
cmap=colormap(colorshere_gray);
c=colorbar
c.Label.String = 'Probability Density (normalized within cluster)';
end

%

if numel(RGB) == 2
for kk = 1:numel(RGB)
hh =figure(480)
%total =RGB{1};
total =(imcomplement(imcomplement(RGB{1})+imcomplement(RGB{2})));%));
%total(total<0) = 0;
%hh=figure(583+kk)
image((total))
%colorbar
axis off
set(hh,'Color','w')
%titlehere = strcat('Merge',strrep(conditionnames{pair(1)},'_','_'),'  ',strrep(conditionnames{pair(2)},'_','_'));
%title(strcat('Merge  ',strrep(conditionnames{pair(1)},'_',''),'  ',strrep(conditionnames{pair(2)},'_','')))
if nargin>4
  %  print('-depsc',strcat(plotdir,filesep,titlehere));
   %     print('-dpng',strcat(plotdir,filesep,titlehere));
end
end
end

end