%% construct a behavioral hierarchy from an input analyssstruct
function [hierarchystruct]= find_sequences_states_demo(analysisstruct,condname_in,params)
annotation_choose = [3:10];
%timescales= [1./16,1./8,1./4,1./2,1,2,4,8,16,32,64]; % timescales in minutes
%timescales= [1./16 1./8 1./4 1./2 1 2 4 8 16];
annotation_choose = [9,10,11,12,22,23,31,44,45,46,29,43];
annotation_choose = find(cellfun(@numel,strfind(analysisstruct.conditionnames,'Baseline1')));
annotation_choose =1;
timescales= [1./16 1./8 1./2  2 5 ];
%timescales =  1./4;

if nargin<3
    params.do_show_pdistmatrix =0;
    params.cases ='real';
    params.decimation_factor = 10;
end

if ~isfield(analysisstruct,'plotdirectory')
    analysisstruct.plotdirectory = 'Y:\Jesse\Data\Plots\Default';
end

%timescales= [ 1./4   2 ];
if ~isfield(params,'timescales')
timescales= [1./8 1./4 1./2 1  2 4 8];
else
    timescales = params.timescales;
end


condnames = condname_in;%{'Baseline1','Baseline2','amph','caff'};%,'PostBi1','PostBi2'
%condnames = {'Baseline1'};
%condnames = {'pre_long_1_JDM25'};%,'ASD'};
%condnames = {'ASD','ASD_control'};
if isfield(params,'corr_threshold')
corr_threshold = params.corr_threshold;
clustercutoff = params.clustercutoff;
else
corr_threshold = 0.3;
clustercutoff = 0.65;
end
chopsize = 20000;
%clustercutoff = 0.65;
%condnames = {'ASD','ephys'};
ratcond_inds = [];
annotation_categories = [];

if ~isfield(params,'doclustering')
params.doclustering = 1;
end


if ~isfield(params,'ratname_sub')
params.ratname_sub = analysisstruct.ratname;
end

for rr = 1:numel(condnames)
    
    ratcond_inds = cat(2,ratcond_inds,find(cellfun(@numel,strfind(analysisstruct.conditionnames,condnames{rr}))));
  
    
    annotation_categories = cat(2,annotation_categories,rr*ones(1,numel(find(cellfun(@numel,strfind(analysisstruct.conditionnames,condnames{rr}))))));
 
end


%% first filter for input ratnames
ratcond_inds = [];
    ratsubset =params.ratname_sub;
annotation_categories = [];
for rr = 1:numel(condnames)

    analysisstruct.conditiontypes_plot = {condnames{rr}};
    analysisstruct = get_matched_conditions(analysisstruct,ratsubset);
ratccondhere = analysisstruct.matchedconds{1};
    ratcond_inds = cat(2,ratcond_inds,ratccondhere);
    annotation_categories = cat(2,annotation_categories,rr*ones(1,numel(ratccondhere)));
end
fprintf('condition names')
analysisstruct.conditionnames(ratcond_inds)

%% specifics for pre ratname
if sum(strcmp(condname_in,'aff'))
ratcond_inds = [];
    ratsubset ={'JDM25','JDM32','JDM33','Vicon8'};
annotation_categories = [];
for rr = 1:numel(condnames)

    analysisstruct.conditiontypes_plot = {condnames{rr}};
    analysisstruct = get_matched_conditions(analysisstruct,ratsubset);
ratccondhere = analysisstruct.matchedconds{1};
    ratcond_inds = cat(2,ratcond_inds,ratccondhere);
    annotation_categories = cat(2,annotation_categories,rr*ones(1,numel(ratccondhere)));
end
end
% %% FOR ASD RESTRICT BY RAT
if strcmp(condname_in,'ASD')
    condnames = {'ASD_control','ASD'};
    
    analysisstruct.conditiontypes_plot = {'ASD_day'};
    ratsubset ={'JDMA1','JDMA2','JDMA3','JDMA4'};
    analysisstruct = get_matched_conditions(analysisstruct,ratsubset);
    matchedcond_ASD = analysisstruct.matchedconds{1};
    
    analysisstruct.conditiontypes_plot = {'ASD_day'};
    ratsubset ={'JDMA5','JDMA6','JDMA8'};
    analysisstruct = get_matched_conditions(analysisstruct,ratsubset);
    matchedcond_ASD2 = analysisstruct.matchedconds{1};
    
    ratcond_inds = cat(2,matchedcond_ASD,matchedcond_ASD2);
    annotation_categories = cat(2,ones(1,numel(matchedcond_ASD)),2*ones(1,numel(matchedcond_ASD2)));
end
%
ratcond_names = analysisstruct.conditionnames(ratcond_inds);
annotation_choose = ratcond_inds;%([1,2,3,7]);

%% -----------------------------------------




%% --------------------------------
nnmf_annotation = cell(1,numel(timescales));
nnmf_score = cell(1,numel(timescales));

annotation_matrix_shifted_conv_agg = cell(1,numel(timescales));
timescalepdist = cell(1,numel(timescales));
seq_inds_cluster_agg= cell(1,numel(timescales));
pdistmatrix  = cell(1,numel(timescales));
outperms = cell(1,numel(timescales));

decimation_factor = params.decimation_factor;
hierarchystruct.decimation_factor=decimation_factor;

%numnnmf_pcs = [5,5];
%1/6 and 25
%% loop over all elements

cases = {'real','markov','shuffled','coarse','neural'};
casehere = find(strcmp(cases,params.cases));

for kkjj=casehere
    
    switch cases{kkjj}
        case 'real'
            base_annotation=[];
            base_annotation_inds=[];
            annotation_indexes=[];
            
            if ~isfield(analysisstruct,'annot_reordered_coarse')                
                for kk =1:numel(annotation_choose)
                    indsuse = 1:numel(analysisstruct.annot_reordered{annotation_choose(kk)});
                    base_annotation = cat(2,base_annotation,reshape(analysisstruct.annot_reordered{annotation_choose(kk)}(indsuse),1,[]));
                    base_annotation_inds = cat(2,base_annotation_inds,annotation_choose(kk)*ones(1,numel(analysisstruct.annot_reordered{annotation_choose(kk)}(indsuse))));
                    annotation_indexes=cat(2,annotation_indexes,find(analysisstruct.condition_inds==annotation_choose(kk))');
                end
            else
                for kk =1:numel(annotation_choose)
                    base_annotation = cat(2,base_annotation,analysisstruct.annot_reordered_coarse{annotation_choose(kk)});
                    base_annotation_inds = cat(2,base_annotation_inds,annotation_choose(kk)*ones(1,numel(analysisstruct.annot_reordered_coarse{annotation_choose(kk)})));
                    annotation_indexes=cat(2,annotation_indexes,find(analysisstruct.condition_inds==annotation_choose(kk))');
                end
            end
            base_annotation_inds_decimation = base_annotation_inds(1:decimation_factor:end);
        case 'markov'
            %% generate markov
            annotation_vector=[];
            base_annotation_inds=[];
            annotation_indexes=[];
            for kk =1:numel(annotation_choose)
                annotation_vector = cat(2,annotation_vector,analysisstruct.annot_reordered{annotation_choose(kk)});
                %% look at all transitions
                indsuse = 1:numel(analysisstruct.annot_reordered{annotation_choose(kk)});
                base_annotation_inds = cat(2,base_annotation_inds,annotation_choose(kk)*...
                    ones(1,numel(analysisstruct.annot_reordered{annotation_choose(kk)}(indsuse))));
                annotation_indexes=cat(2,annotation_indexes,find(analysisstruct.condition_inds==annotation_choose(kk))');
            end
            maxval =max(unique(annotation_vector)); %include null/rest as 0
            transprob = zeros(maxval+1,maxval+1);
            statefreq = arrayfun(@(x) numel(find(annotation_vector == x)),0:maxval);
            
            
            for mm = 1:numel(annotation_vector)-1
                transprob(annotation_vector(mm)+1,annotation_vector(mm+1)+1) = transprob(annotation_vector(mm)+1,annotation_vector(mm+1)+1)+...
                    1./statefreq(annotation_vector(mm)+1);
            end
            %get rid of the diagonal
            transprob_offdiag = transprob-diag(diag(transprob));
            %normalize by the size of the diagonal
            exit_prob = bsxfun(@rdivide,transprob_offdiag,(1-diag(transprob)));
            mc = dtmc(exit_prob);
            try
            base_annotation = simulate(mc,numel(annotation_vector));
            catch ME
                            base_annotation = simulate(mc,numel(annotation_vector));
            end
            % other housekeeping indicies
            base_annotation = reshape(base_annotation,1,[]);
            base_annotation_inds_decimation = base_annotation_inds(1:decimation_factor:end);
            
        case 'coarse'
            base_annotation=[];
            base_annotation_inds=[];
            for kk =1:numel(annotation_choose)
                base_annotation = cat(2,base_annotation,analysisstruct.full_annotation_reduced{annotation_choose(kk)});
                base_annotation_inds= cat(2,base_annotation_inds,annotation_choose(kk)*ones(1,numel(analysisstruct.full_annotation_reduced{annotation_choose(kk)})));
            end
            base_annotation_inds_decimation = base_annotation_inds(1:decimation_factor:end);
            %% clustering behavior by neural data
         case 'neural'
            base_annotation=[];
            base_annotation_inds=[];
            annotation_indexes=[];
                  
                for kk =1:numel(annotation_choose)
                    
                    base_annotation = cat(2,base_annotation,params.neural_matrix{annotation_choose(kk)});
                    base_annotation_inds = cat(2,base_annotation_inds,annotation_choose(kk)*ones(1,size(params.neural_matrix{annotation_choose(kk)},1)));
                    annotation_indexes=cat(2,annotation_indexes,find(analysisstruct.condition_inds==annotation_choose(kk))');
                end           
            base_annotation_inds_decimation = base_annotation_inds(1:decimation_factor:end);
    end
end

%% get similarity matrix over timescales
sort_dendrogram = 0;
do_save_videos = 0;



%% matricies for saving
clustered_behavior = cell(1,numel(timescales));
clustered_sequences = cell(1,numel(timescales));
sequences_inds = cell(1,numel(timescales));

hierarchystruct = struct();
hierarchystruct.annotation_categories = annotation_categories;
%hierarchy_struct.conditionnames
%hierarchy_struct.conditionnumbers
%hierarchystruct.behavior_annotations

for ll = 1:numel(timescales)
    circshiftval = round((300./(analysisstruct.tsnegranularity))*60*timescales(ll));
    %if (ll==1)395
    annotation_use = base_annotation;
    if ~strcmp(cases{casehere},'neural')
    annotation_matrix_shifted = zeros(max((base_annotation)),numel(base_annotation));
    for rr = reshape((unique(base_annotation(base_annotation>0))),1,[])
        annotation_matrix_shifted(rr,find(annotation_use==rr)) = 1;
    end
    annotation_matrix_shifted(find(sum(annotation_matrix_shifted,2)==0),:) =[];
    
    % get PDF of behaviors at each timepoint
    annotation_matrix_shifted_conv = conv2(annotation_matrix_shifted,ones(1,circshiftval),'same');
    annotation_matrix_shifted_conv_agg = bsxfun(@rdivide,annotation_matrix_shifted_conv,nansum(annotation_matrix_shifted_conv,1));
    annotation_matrix_shifted_conv_agg(isnan(annotation_matrix_shifted_conv_agg)) = 0;
    else
        base_sz = size(base_annotation);
        if (base_sz(2)<base_sz(1))
        base_annotation = base_annotation';
        end
       annotation_matrix_shifted_conv_agg = conv2(base_annotation,ones(1,circshiftval),'same');
    end
    do_show_pdistmatrix =  params.do_show_pdistmatrix ;
    if do_show_pdistmatrix
        timescalepdist = 1-pdist(annotation_matrix_shifted_conv_agg(:,1:decimation_factor:end)','correlation');
        valplot_nonbiased = squareform(timescalepdist);
        figure(434)
        subplot(1,numel(timescales),ll)
        imagesc(valplot_nonbiased)
        colormap(cat(1,ones(3,3),othercolor('PuRd7',25)))
        
        valplot_nonbiased(find(tril(ones(size(valplot_nonbiased)),round(2*(circshiftval./decimation_factor))))) = nan;
        
        clear valplot_nonbiased
        clear timescalepdist
        
    end

    %% begin breaking things up
    %
    %  circshiftval = round((300./(analysisstruct.tsnegranularity))*60*timescales(ll));
    pixelnums = round(2*(circshiftval./decimation_factor));
    size_threshold = pixelnums;
    annotation_matrix_shifted_conv_agg = annotation_matrix_shifted_conv_agg(:,1:decimation_factor:end);
    %% break up watershed to get over memory complexity issue
    numchops = ceil(size( annotation_matrix_shifted_conv_agg,2)./chopsize);
    density_cc = struct('NumObjects',[0],'PixelIdxList',{cell(0,0)});
    limhere = size(annotation_matrix_shifted_conv_agg,2);
    corr_hist_range = -1:0.025:1;
    histsum = zeros(1,numel(corr_hist_range));
    histsum_crosscomp = zeros(numel(annotation_choose),numel(annotation_choose),numel(corr_hist_range));
    
    for jj=1:numchops
        for kk=1:numchops
            fprintf('chopping up distance matrix for %f and %f of %f \n',jj,kk,numchops);
            vals_1 = bsxfun(@plus,1:chopsize,(jj-1)*chopsize);
            vals_1(vals_1>limhere) = [];
            vals_2 = bsxfun(@plus,1:chopsize,(kk-1)*chopsize);
            vals_2(vals_2>limhere) = [];
            pdistmatrixhere= 1-pdist2(annotation_matrix_shifted_conv_agg(:,vals_1)',...
                annotation_matrix_shifted_conv_agg(:,vals_2)','correlation');
            %% STILL HAVE TO NAN OUT THE DIAGONAL, if on the diagonal
            % valplot_nonbiased=1-pdist(vals_1*decimation_factor);
            if jj==kk
                pdistmatrixhere(find(tril(ones(size(pdistmatrixhere)),round(2*(circshiftval./decimation_factor))))) = nan;
            end
            
            %store values for later on
            [nn,xx] = hist(pdistmatrixhere(~isnan(pdistmatrixhere)),corr_hist_range);
            histsum = histsum+nn;
            
            %% also do for indiv conditions
            for kjk =1:numel(annotation_choose)
                for jjkk =1:numel(annotation_choose)
                    ind1 = find(ismember(vals_1,find(base_annotation_inds_decimation == annotation_choose(kjk))));
                    ind2 = find(ismember(vals_2,find(base_annotation_inds_decimation == annotation_choose(jjkk))));
                    if numel(ind1) && numel(ind2)
                        vals_examine_here = pdistmatrixhere(ind1,ind2);
                        [histsum_crosscomp(kjk,jjkk,:),~] = hist(vals_examine_here(~isnan(vals_examine_here)),corr_hist_range);
                    end
                end
            end
            
            % threshold for making clusters
            if numel(corr_threshold)==1
            pdistmatrixhere(pdistmatrixhere<corr_threshold) = 0;
            else
                  pdistmatrixhere(pdistmatrixhere<corr_threshold(ll)) = 0;          
            end
            pdistmatrixhere(isnan(pdistmatrixhere)) = 0;
            
            L = watershed(pdistmatrixhere);%valplot_nonbiased(vals_1,vals_2));
            density_cc_temp = bwconncomp(L);
            numinds = cellfun(@numel,density_cc_temp.PixelIdxList);
            [~,maxind] = max(numinds);
            goodobj =0;
            goodinds = setxor(find(numinds>=size_threshold),maxind);
            for rr = goodinds
                %     if rr~=maxind && numel( density_cc_temp.PixelIdxList{rr})>size_threshold
                [sub1,sub2] = ind2sub(size(L),density_cc_temp.PixelIdxList{rr});
                density_cc.PixelIdxList{end+1} = sub2ind([size(annotation_matrix_shifted_conv_agg,2) size(annotation_matrix_shifted_conv_agg,2)],...
                    bsxfun(@plus,sub1,(jj-1)*chopsize),bsxfun(@plus,sub2,(kk-1)*chopsize));
                goodobj = goodobj+1;
                %    end
            end
            
            
            % L(find(tril(ones(size(valplot)),pixelnums))) = 0;
            density_cc.NumObjects = density_cc.NumObjects+goodobj;
        end
    end
    clear vals_examine_here pdistmatrixhere density_cc_temp
    
    figure(5454)
    subplot(1,numel(timescales),ll)
    bar(corr_hist_range,histsum./sum(histsum))
    box off
    ylim([0 0.2])
    %set(gca,'YScale','log')
    %
    
    if numel(annotation_choose)>=2
        colorshere = othercolor('Mrainbow',numel(unique(annotation_categories)));
        figure(54588)
        subplot(1,numel(timescales),ll)
        for kjk =1:numel(annotation_choose)
            %   plot(corr_hist_range,squeeze(histsum_crosscomp(kjk,kjk,:))./sum(squeeze(histsum_crosscomp(kjk,kjk,:))))
            %    hold on
            %    box off
            %   ylim([0 0.2])
        end
        for nnk = unique(annotation_categories)
            av_val = zeros(numel(corr_hist_range),1);
            for kjk = find(annotation_categories==nnk)
                av_val = av_val+squeeze(histsum_crosscomp(kjk,kjk,:))./sum(squeeze(histsum_crosscomp(kjk,kjk,:)));
            end
            plot(corr_hist_range,av_val./numel(find(annotation_categories==nnk)),'linewidth',2,'Color',colorshere(nnk,:))
            hold on
        end
        legend(condnames) %cat(2,{'baseline'},{'asd'}))
        
        %plot after to improve the legend
        for nnk = unique(annotation_categories)
            for kjk = find(annotation_categories==nnk)
                plot(corr_hist_range,squeeze(histsum_crosscomp(kjk,kjk,:))./sum(squeeze(histsum_crosscomp(kjk,kjk,:))),'Color',colorshere(nnk,:));
            end
        end
        hold off
        %    if ll==2
        %  legend(cat(2,analysisstruct.conditionnames(annotation_choose),{'baseline'},{'asd'}))
        
        % end
    end
    
    if numel(annotation_choose)>=2
        colorshere = othercolor('Mrainbow',numel(unique(annotation_categories)));
        figure(54589)
        subplot(1,numel(timescales),ll)
        for nnk = unique(annotation_categories)
            av_val = zeros(numel(corr_hist_range),1);
            for kjk = find(annotation_categories==nnk)
                av_val = av_val+squeeze(histsum_crosscomp(kjk,kjk,:))./sum(squeeze(histsum_crosscomp(kjk,kjk,:)));
            end
            semilogy(corr_hist_range,av_val./numel(find(annotation_categories==nnk)),'linewidth',2,'Color',colorshere(nnk,:))
            hold on
        end
        legend(condnames)
        box off
        legend boxoff
        ylabel('fraction of all frames')
        ylim([0.0005 1])
        hold off
    end
    
    %
    %      if ll==numel(timescales)
    %    print('-depsc',strcat(analysisstruct.plotdirectory,filesep,'Sequence_usage.epsc'))
    %        print('-dpng',strcat(analysisstruct.plotdirectory,filesep,'Sequence_usage.png'))
    %     end
    
    
    %% get the integrated correlation over this timescale
    integrated_correlation_comparison = zeros(numel(annotation_choose),numel(annotation_choose));
    for kjk =1:numel(annotation_choose)
        for lkl =1:numel(annotation_choose)
            integrated_correlation_comparison(kjk,lkl) = ...
                0.5*( nansum(corr_hist_range'.*squeeze(histsum_crosscomp(kjk,lkl,:))./nansum(squeeze(histsum_crosscomp(kjk,lkl,:))))+...
                nansum(corr_hist_range'.*squeeze(histsum_crosscomp(lkl,kjk,:))./nansum(squeeze(histsum_crosscomp(lkl,kjk,:)))));
        end
    end
    
    figure(55)
    subplot(1,numel(timescales),ll)
    imagesc(integrated_correlation_comparison)
    set(gca,'XTick',1:numel(annotation_choose),'XTickLabel',analysisstruct.conditionnames(annotation_choose))
    xtickangle(90)
    set(gca,'YTick',1:numel(annotation_choose),'YTickLabel',analysisstruct.conditionnames(annotation_choose))
    
    %% also get the averages cross for all conditions
    crossvals = zeros(numel(unique(annotation_categories)),numel(unique(annotation_categories)));
    crossvals_std = zeros(numel(unique(annotation_categories)),numel(unique(annotation_categories)));
    nameshere = cell(numel(unique(annotation_categories)),numel(unique(annotation_categories)));
    for nnk = unique(annotation_categories)
        for nnj = unique(annotation_categories)
            crossvals(nnk,nnj) = nanmean(nanmean(integrated_correlation_comparison(find(annotation_categories==nnk),find(annotation_categories==nnj))));
            valshere = integrated_correlation_comparison(find(annotation_categories==nnk),find(annotation_categories==nnj));
            crossvals_std(nnk,nnj) = nanstd(valshere(:))./sqrt(numel(valshere));
            nameshere{nnk,nnj} = strcat(condnames{nnk},{' '},condnames{nnj});
            nameshere{nnk,nnj} =     nameshere{nnk,nnj}{1}
        end
    end
    figure(56)
    subplot(1,numel(timescales),ll)
    
    indsuse = find(triu(ones(size(crossvals))));
    bar(crossvals(indsuse))
    hold on
    errorbar(1:numel(indsuse),crossvals(indsuse),crossvals_std(indsuse),'k','Linestyle','none')
    set(gca,'XTick',1:numel(indsuse),'XTickLabels',nameshere(indsuse))
    ylabel('Correlation between conditions')
    xtickangle(90)
    
    %% save to the hierarchy structure
    hierarchystruct.correlation_histograms{ll} = histsum;
    hierarchystruct.histxrange{ll} = corr_hist_range;
    hierarchystruct.correlation_comp_histograms{ll} = histsum_crosscomp;
    hierarchystruct.integrated_correlation_histograms{ll} =integrated_correlation_comparison;
    
    %% get the sclusters for these
    numinds = cellfun(@numel,density_cc.PixelIdxList);
    [~,indbackground] = max(numinds);
    bad_inds = cat(2,indbackground,find(numinds <= size_threshold));
    
    
    good_inds = setxor(1:density_cc.NumObjects,bad_inds);
    mean_seq_val = zeros(size(annotation_matrix_shifted_conv_agg,1),numel(good_inds));
    seq_inds = cell(numel(good_inds),2);
    
    for lkl=1:numel(good_inds)
        [seq_inds{lkl,1},seq_inds{lkl,2}] = ind2sub([size(annotation_matrix_shifted_conv_agg,2) size(annotation_matrix_shifted_conv_agg,2)],...
            density_cc.PixelIdxList{good_inds(lkl)});
        mean_seq_val(:,lkl) = mean(annotation_matrix_shifted_conv_agg(:,[seq_inds{lkl,1} seq_inds{lkl,2}]),2);
    end
    sequences_inds{ll} = seq_inds;
    clustered_sequences{ll} = mean_seq_val;
end

if casehere ==2
    fprintf('doing markov only \n')
    %return;
end
fprintf('now starting clustering \n')
equivclass_fraction_agg = cell(1,numel(timescales));
normshades_agg = cell(1,numel(timescales));
behshades_agg = cell(1,numel(timescales));
behavior_composition_vector = cell(1,numel(timescales));

if strcmp(cases{casehere},'neural')
                             base_annotation = reshape(analysisstruct.annot_reordered{annotation_choose(kk)},1,[]);
end

if ~params.doclustering
    return
end
for ll = 1:numel(timescales)
    %% clustering step
    fprintf('correlation based clustering over %f points \n',size(clustered_sequences{ll},2))
    
    if ~strcmp(cases{casehere},'neural')
    behsuse = find(nanmean(clustered_sequences{ll},2)>0.0001);
    else
        behsuse = 1:size(clustered_sequences{ll},1);
    end
    
    if numel(clustercutoff)>1
        clustercutoff_here = clustercutoff(ll);
    else
        clustercutoff_here = clustercutoff;
    end
    
    nsample = 55000;
    num_obs = size(clustered_sequences{ll}(behsuse,:),2);
    tic
    fprintf('Number of clusters %f doing knn reembed or not \n',num_obs);
    if num_obs<nsample
                fprintf('Clustering (no reembed needed) \n')
            Z = linkage(clustered_sequences{ll}(behsuse,:)','average','correlation');
                T = cluster(Z,'cutoff',clustercutoff_here,'criterion','distance','depth',3);
    else
        
        nclust = 300;
        fprintf('importance sampling \n')
        iddx = kmeans(clustered_sequences{ll}(behsuse,:)',nclust);
        [un_el,~,ic] = unique(iddx);
                            counts = accumarray(ic,1);
n_guide_points = [];
n_subsample = floor(nsample./nclust);
for nn = unique(iddx)'
            n_guide_points = cat(1,n_guide_points,randsample(find(iddx==nn),min(counts(nn),n_subsample)));
end
                            
       % n_guide_points = randsample(1:num_obs,nsample);
        n_oos = setxor(1:num_obs,n_guide_points);
        
        fprintf('Clustering \n')
         Z = linkage(clustered_sequences{ll}(behsuse,n_guide_points)','average','correlation');
         T_is = cluster(Z,'cutoff',clustercutoff_here,'criterion','distance','depth',3);
        fprintf('Knn rembed \n')
         idx = knnsearch(clustered_sequences{ll}(behsuse,n_guide_points)',...
             clustered_sequences{ll}(behsuse,n_oos)');
         T_oos = T_is(idx);
         
         T= zeros(1,num_obs);
         T(n_guide_points) = T_is;
         T(n_oos) = T_oos;
    end
    toc
 
   %  clustercutoff=0.4
   % T = cluster(Z,'cutoff',clustercutoff,'criterion','distance','depth',3);
    color = Z(end-4,3)-eps;
      

     %dendrogram(Z,25)
     axis off

    % this wasnt great
    %T = clusterdata(clustered_sequences{ll}','criterion','distance','Cutoff',1.9,'linkage','ward','Depth',3,'distance','euclidean','SaveMemory','on');
    %figure(44)
    % imagesc(squareform(pdist(mean_seq_val','euclidean')))
    %% setup the loop over clusters
    seq_inds_cluster_agg{ll} = cell(1,numel(unique(T)));
    conditions_using_seq =cell(1,numel(unique(T)));
    base_annot_clustered = zeros(1,numel(base_annotation));
    base_annot_vectors = zeros(numel(unique(annotation_choose)),numel(unique(T)));
   
    %% get the number of 'good' clusters
    [un_el,~,ic] = unique(T);
    counts_seq = accumarray(ic,1);
    goodclust = find(counts_seq>20);
    fprintf('numel clust (min 20 fr) %f for timescale %f \n',  numel(goodclust),timescales(ll)*60);
    [~,gcsorted] = sort(counts_seq(goodclust),'DESCEND');
    
    if ll==1
    num_ex = 5;
    else
        num_ex=10;
    end
    indsgood=[];
    for lk = goodclust'
    indsgood = cat(1,indsgood,randsample(find(ismember(T,lk)),num_ex,'true'));
    end
     
     figure(788+ll)
%[H,T,perm] = dendrogram(Z, 25, 'colorthreshold', color);
    Ztwo = linkage(clustered_sequences{ll}(behsuse,indsgood)','average','correlation');

[H,~,~] = dendrogram(Ztwo, 'colorthreshold', 0.7);
    
    %  sum_corr = sum(valplot_nonbiased,2);
    % unstruct_frames = find(sum_corr<30);
    %      unique_total_inds = unique(round(reshape(bsxfun(@plus,unstruct_frames*decimation_factor,-decimation_factor./2:1:decimation_factor./2),[],1)));
    %  unique_total_inds(unique_total_inds>numel(base_annotation)) = numel(base_annotation);
    %  M=animate_markers_aligned_fullmovie(analysisstruct.mocapstruct_reduced_agg{1},find(base_annot_clustered==6))
    
    allbehshere = reshape(unique(base_annotation),1,[]);
    do_save_videos=0;
    beh_shades_agg = zeros(numel(gcsorted),analysisstruct.density_objects+1);
    behs_in_cluster = [];
    for clusterhere = 1:numel(gcsorted)
        clusterinds = find(T==goodclust(gcsorted(clusterhere)));
        
        %% seq inds
        % ind_cell = cell(0,0);
        % for kk = clusterinds'
        %    seq_inds_cluster =  cat(1,seq_inds{(kk),1},seq_inds{kk,2});
        %       ind_cell{end+1}=unique(reshape(bsxfun(@plus,seq_inds_cluster*decimation_factor,-decimation_factor./2:1:decimation_factor./2),[],1));
        % end
        seq_inds_cluster = unique(cat(1,sequences_inds{ll}{clusterinds,1},...
            sequences_inds{ll}{clusterinds,2}));
        max(seq_inds_cluster)
        unique_total_inds = unique(round(reshape(bsxfun(@plus,...
            seq_inds_cluster*decimation_factor,-decimation_factor./2:1:decimation_factor./2),[],1)));
        unique_total_inds(unique_total_inds>numel(base_annotation_inds)) = numel(base_annotation_inds);
        conditions_using_seq{clusterhere} = unique(base_annotation_inds(unique_total_inds));
        base_annot_clustered(unique_total_inds) = clusterhere;
        for kz = unique(annotation_choose)
            base_annot_vectors(find(annotation_choose == kz),clusterhere) = ...
                numel(intersect( unique_total_inds,find(base_annotation_inds == kz)))./numel(find(base_annotation_inds == kz));
        end
        [fulloutput,indivbouts] = fillannotationgaps(unique_total_inds',1,0);
        seq_inds_cluster_agg{ll}{clusterhere} = seq_inds_cluster;
        %% plot cluster colored tsne
        beh_shades = zeros(1,analysisstruct.density_objects+1);
        
        [un_el,~,ic] = unique(base_annotation(unique_total_inds));
        counts = accumarray(ic,1);
        beh_shades(un_el+1) = counts./sum(counts);
        beh_shades(1) = [];
        beh_shades_agg(clusterhere,1:numel(beh_shades)) =beh_shades;
        %  beh_shades(beh_shades>0.25) = 0.25;
        %       beh_shades(beh_shades<0.01) = 0.0;
        
        normshades = log2(beh_shades);
        normshades(isinf(normshades)) = -13;
        normshades_agg{ll}{clusterhere} = normshades;
        behshades_agg{clusterhere} = beh_shades_agg;
        
        
        
        figure(395+ll)
        h1=subplot_tight(7,3,clusterhere)
        % plot_shaded_tsne(analysisstruct,normshades,[],1,h1)
        %         plot_shaded_tsne(analysisstruct,normshades,[],1,h1)
        plot_shaded_tsne_absolute_alpha(analysisstruct,normshades,[],1,h1,1)
        axis equal
        axis off
        colorbar off
        %  title(num2str(counts_seq(goodclust(gcsorted(clusterhere)))))
        caxis([0.01 0.25])
        if strcmp(cases{kkjj},'coarse')
            fprintf('for cluster %f frac move: %f frac fast %f frac rest %f \n',clusterhere, beh_shades(end-2),beh_shades(end-1),beh_shades(end));
        end
        if clusterhere == numel(gcsorted)
            mkdir(analysisstruct.plotdirectory)
            set(gcf,'renderer','opengl')
            print('-depsc',strcat(analysisstruct.plotdirectory,filesep,'sequence_maps',num2str(ll),'.epsc'))
            print('-dpng',strcat(analysisstruct.plotdirectory,filesep,'sequence_maps',num2str(ll),'.png'))
        end
        
        
        analysisstruct.clusternames_sorted=analysisstruct.clusternames;
        
        figure(495+ll)
        h1=subplot_tight(7,3,clusterhere)
        uniquevals = unique(analysisstruct.clusternames_sorted );
        unique_equivclass = cell(1,numel(uniquevals));
        equivclass_fraction = zeros(1,numel(uniquevals));
        for zz = 1:numel(uniquevals)
            unique_equivclass{zz} = find(strcmp(analysisstruct.clusternames_sorted ,uniquevals{zz}));
            equivclass_fraction(zz) = sum(squeeze(beh_shades_agg(clusterhere,(unique_equivclass{zz}))));
        end
        bar(equivclass_fraction(analysisstruct.behavior_order_numbers),'FaceColor',[0.5 .5 .5],'EdgeColor','none')
        box off
        if clusterhere>(numel(gcsorted)-2)
            set(gca,'XTick',1:numel(equivclass_fraction),'XTickLabels',uniquevals(analysisstruct.behavior_order_numbers))
            xtickangle(90)
            
            figure(499)
            bar(equivclass_fraction(analysisstruct.behavior_order_numbers),'FaceColor',[0.5 .5 .5],'EdgeColor','none')
            set(gca,'XTick',1:numel(equivclass_fraction),'XTickLabels',uniquevals(analysisstruct.behavior_order_numbers))
            xtickangle(90)
            box off
        end
        equivclass_fraction_agg{ll}{clusterhere} = equivclass_fraction;
        hierarchystruct.uniquenames = uniquevals;
        %goodclust(gcsorted(clusterhere))
        
        %  animate_markers_aligned_fullmovie(analysisstruct.mocapstruct_reduced_agg{1},...
        %   round(unique(reshape(bsxfun(@plus,seq_inds_cluster*decimation_factor,-decimation_factor*2:1:decimation_factor*2),[],1))))
        if do_save_videos
            figure(388)
            clf
            plotfolder = strcat('Y:\Jesse\Data\Motionanalysis_captures\Cluster_videos\20181206_timescaleplots_',num2str(timescales(ll)),'\');
            mkdir(plotfolder)
            [~, biggestindiv] = sort(cellfun(@numel,indivbouts),'DESCEND')
            M=plot_multi_clusters(analysisstruct.mocapstruct_reduced_agg{1},indivbouts(biggestindiv(1:6)),numel(indivbouts{1}))
            v = VideoWriter(strcat(plotfolder,filesep,'clusteredmovie_',num2str((clusterhere))),'MPEG-4');
            open(v)
            writeVideo(v, M)
            close(v)
        end
        
        %% get behaviors in cluster
        %% get the number of indiv bouts for each behavior
        
        boutlist = [];
        for rr = 1:numel(indivbouts)
            boutlist = cat(1,boutlist,rr*ones(numel(indivbouts{rr}),1));
        end
        boutnums = zeros(1,numel(allbehshere));
        for llk = allbehshere
            boutnums(find(allbehshere==llk)) = numel(unique(boutlist(find(base_annotation(cat(2,indivbouts{:}))==llk))));
        end
        
        %behaviors in at least 1/8th
        behs_in_cluster{clusterhere} = allbehshere(find(boutnums>=4));%floor(numel(indivbouts)./8)));
        
    end
    
    %
    figure(501+ll)
    cmaphere = cat(1,ones(3,3),othercolor('PuRd9',50));
    equivclass_fraction = zeros(size(beh_shades_agg,1),numel(uniquevals));
    for zz = 1:numel(uniquevals)
        %  unique_equivclass{zz} = find(strcmp(analysisstruct.clusternames_sorted ,uniquevals{zz}));
        equivclass_fraction(:,zz) = sum(beh_shades_agg(:,unique_equivclass{zz} ),2);
    end
    imagesc(equivclass_fraction')
    colormap(cmaphere)
    colorbar
    set(gca,'YTick',1:numel(equivclass_fraction),'YTickLabels',uniquevals)
    xlabel('State ID')
    box off
    title('State Composition')
    clustered_behavior{ll} = base_annot_clustered;
    
    %% compare the behavioral usage across conditions
    base_condition = 1;
    % for jkj=1:3
    behavior_vector = zeros(numel(annotation_choose),max(clustered_behavior{ll})+1);
    for kk =1:numel(annotation_choose)
        beh_use_vec =  clustered_behavior{ll}(find(base_annotation_inds ==annotation_choose(kk)));
        [un_el,~,ic] = unique(beh_use_vec);
        counts = accumarray(ic,1)./numel(beh_use_vec);
        %% get the finger print for each type
        behavior_vector(kk,un_el+1) = counts;
        fprintf('amount of time unstructured %s is %f \n',analysisstruct.conditionnames{annotation_choose(kk)},behavior_vector(kk,1));
        
    end
    %get average behavior vectors
    behvector_av = zeros(numel(unique(annotation_categories)),size(behavior_vector,2));
    for nnk = unique(annotation_categories)
        behvector_av(nnk,:) = nanmean(behavior_vector(find(annotation_categories==nnk),:),1);
    end
    
    figure(333)
    subplot(1,numel(timescales),ll)
    bar(behavior_vector(:,1))
    set(gca,'XTick',1:numel(annotation_choose),'XTickLabel',analysisstruct.conditionnames(annotation_choose))
    ylabel('percent of time unstructured')
    xtickangle(90)
    
    figure(324)
    subplot(1,numel(timescales),ll)
    imagesc(behavior_vector(:,2:end))
    set(gca,'YTick',1:numel(annotation_choose),'YTickLabel',strrep(analysisstruct.conditionnames(annotation_choose),'_',' '))
    title('behavioral usage vector')
    colormap(cat(1,ones(3,3),othercolor('PuRd7',25)))
    if numel(behavior_vector)>1
    caxis([0 0.9*max(max(behavior_vector(:,2:end)))])
    end
    % set(gca,'FontSize',18)
    behavior_composition_vector{ll} = behavior_vector;
    xlabel('Behavioral States')
    box off
    if ll==numel(timescales)
        print('-depsc',strcat(analysisstruct.plotdirectory,filesep,'Sequence_usage_overrats.epsc'))
        print('-dpng',strcat(analysisstruct.plotdirectory,filesep,'Sequence_usage_overrats.png'))
    end
    
    
    
    figure(3245)
    subplot(1,numel(timescales),ll)
    thresh = [0.01 0.01 0.01 0.01 0.01 0.01 0.01];
    goodinds = find(sum(behvector_av,1)>thresh(ll));
    goodinds = setxor(goodinds,intersect(goodinds,0));
    imagesc(behvector_av(:,goodinds))
    colormap(cat(1,ones(3,3),othercolor('PuRd7',25)))
    box off
    caxis([0 0.9*max(max(behvector_av(:,goodinds)))])
    xlabel('Behavioral States')
    c = colorbar
    c.Label.String = 'State Probability Density';
    c.Label.FontSize = 18;
    c.FontSize = 16;
    set(gca,'FontSize',16)
    set(gca,'YTick',1:size(behvector_av(:,:),1))
    set(gcf,'Renderer','Painters')
    print('-depsc',strcat(analysisstruct.plotdirectory,filesep,'behavioral_usage_heatmap_USETHIS.epsc'))
    print('-dpng',strcat(analysisstruct.plotdirectory,filesep,'behavioral_usage_heatmap_USETHIS.png'))
    hierarchystruct.goodstateinds{ll} = goodinds;
    
    figure(327)
    colormap(cat(1,ones(3,3),othercolor('PuRd7',25)))
    caxis([0 0.25]) %0.9*max(max(behavior_vector(:,2:end)))
    axis off
    c = colorbar
    c.Label.String = 'State Probability Density';
    c.Label.FontSize = 18;
    c.FontSize = 16;
    if ll==numel(timescales)
        print('-depsc',strcat(analysisstruct.plotdirectory,filesep,'behavioral_density_axis_USETHIS.epsc'))
        print('-dpng',strcat(analysisstruct.plotdirectory,filesep,'behavioral_density_axis_USETHIS.png'))
    end
    
    %% look across conditions to see if any difference across conditions
    behvals = zeros(numel(unique(annotation_categories)),size(behavior_vector,2));
    behstd = zeros(numel(unique(annotation_categories)),size(behavior_vector,2));
    
    for nnk = unique(annotation_categories)
        behvals(nnk,:) = nanmean(behavior_vector(find(annotation_categories==nnk),:),1);
        behstd(nnk,:) = nanstd(behavior_vector(find(annotation_categories==nnk),:),[],1);
    end
    
    %
    if size(behvals,1)>1
        statescomare = [2,1]
        goodbeh = find((sum(behvals(statescomare,:),1)>0.05));
        goodbeh(1) = [];
        fractionalchange = (behvals(statescomare(1),goodbeh)-behvals(statescomare(2),goodbeh))./behvals(statescomare(2),goodbeh);
        figure(345)
        subplot(1,numel(timescales),ll)
        bar(fractionalchange)
        hold on
        errorbar(1:numel(goodbeh),fractionalchange,...
            nanmean(behstd(statescomare,goodbeh),1)./(sqrt(14)*behvals(statescomare(2),goodbeh)),'k','Linestyle','none')
        hold off
        box off
        set(gca,'FontSize',18)
        xlabel('Behavioral State')
        ylabel('Fractional Change in State Frequency')
        if ll==numel(timescales)
            print('-depsc',strcat(analysisstruct.plotdirectory,filesep,'Sequence_usage.epsc'))
            print('-dpng',strcat(analysisstruct.plotdirectory,filesep,'Sequence_usage.png'))
        end
        
        
        [sortvals,sortedinds] = sort((behvals(2,:)-behvals(1,:)));
        sortedinds(sortedinds==1) = [];
        figure(33)
        subplot(1,numel(timescales),ll)
        plot(behvals(:,sortedinds)','linewidth',2)
        legend(condnames)
        box off
        legend boxcloff
        xlabel('sequence ID')
        ylabel('State frequency')
        title('Behavioral State Usage')
        if ll==numel(timescales)
            print('-depsc',strcat(analysisstruct.plotdirectory,filesep,'Sequence_usage.epsc'))
            print('-dpng',strcat(analysisstruct.plotdirectory,filesep,'Sequence_usage.png'))
        end
        
        figure(3244)
        subplot(1,numel(timescales),ll)
        imagesc(behavior_vector(:,sortedinds))
        set(gca,'YTick',1:numel(annotation_choose),'YTickLabel',analysisstruct.conditionnames(annotation_choose))
        title('behavioral usage vector (sorted)')
        
        
        figure(33333)
        subplot(numel(timescales),1,ll)
        imagesc(equivclass_fraction(sortedinds-1,:)')
        colormap(cmaphere)
        colorbar
        set(gca,'YTick',1:numel(equivclass_fraction),'YTickLabels',uniquevals)
        xlabel('State ID')
        box off
        title('State Composition')
        %find the name
        %ratcond_names{14}
        condex = ratcond_inds(end-2);
        seq_ex =   sortedinds(min(7,numel(sortedinds)))-1;%       sortedinds(4);
        vals_seq = zeros(1,numel(ratcond_inds));
        for lkk = 1:numel(ratcond_inds)
            framesuse = find(ismember(find(base_annotation_inds==ratcond_inds(lkk)),find(clustered_behavior{ll}==seq_ex)));
            vals_seq(lkk) = numel(framesuse);
        end
        [~,maxind] = max(vals_seq);
        framesuse = find(ismember(find(base_annotation_inds==ratcond_inds(maxind)),find(clustered_behavior{ll}==seq_ex)));
        
        %
        %    M=animate_markers_aligned_fullmovie(analysisstruct.mocapstruct_reduced_agg{ratcond_inds(maxind)},...
        %       framesuse(1:min(1000,numel(framesuse)) ))
        %
        
        %       v = VideoWriter(strcat(analysisstruct.plotdirectory,filesep,'clusteredmovie_asd.avi'),'MPEG-4');
        %         open(v)
        %         writeVideo(v, M)
        %         close(v)
        %
        
        %Dwell time in states
        
        dwell_times = cell(numel(unique(annotation_categories)),1);%numel(analysisstruct.sorted_clust_ind));
        dwell_times_av = zeros(numel(unique(annotation_categories)),size(behavior_vector,2));
        for lk = unique(annotation_categories)
            dwell_times{lk} =zeros(numel(find(annotation_categories==lk)),size(behavior_vector,2));
            subsethere = find(annotation_categories==lk);
            for mm = 1:numel(subsethere)
                for ms =0:size(behavior_vector,2)-1
                    indexhere =     find(base_annotation_inds==ratcond_inds(subsethere(mm)));
                    aa = zeros(1,numel(clustered_behavior{ll}(indexhere)));
                    aa(clustered_behavior{ll}(indexhere) == ms) = 1;
                    rr = bwconncomp(aa);
                    dwell_times{lk}(mm,ms+1) = mean(cellfun(@numel,rr.PixelIdxList));
                end
            end
            dwell_times_av(lk,:) = nanmean(  dwell_times{lk},1);
        end
        
        figure(3404)
        plot(dwell_times_av(:,sortedinds)')
        
    end
    
    
    
    
    %% do the behavior across all rats to look at individuality
    behvals = zeros(numel(unique(base_annotation_inds)),size(behavior_vector,2)-1);
    
    for nnk = 1:numel(unique(base_annotation_inds))
        behvals(nnk,:) = nanmean(behavior_vector(nnk,2:end),1);
    end
    
    
    figure(31)
    subplot(1,numel(timescales),ll)
    plot(behvals')
    legend(ratcond_names)
    xlabel('sequence ID')
    ylabel('Behavioral frequency')
    title('Behavioral Usage Across All individual Conditions')
    
    
    %% look at the number of conditions between behaviors
    use_thresh = 0.01;
    behval_thresh = behvals;
    behval_thresh(behval_thresh<use_thresh) = 0;
    behval_thresh(behval_thresh>0) = 1;
    num_beh_ea_cond = sum(behval_thresh,1);
    figure(30)
    [n,x] = ecdf(num_beh_ea_cond);
    plot(x,n)
    title('Sequence use across conditions')
    xlabel('Number of conditions')
    ylabel('Fraction of sequences')
    
    % similarity of conditions
    %% correlation of behavioral usage vectors across conditions
    condition_sequence_similarity = pdist(behavior_vector(:,2:end),'correlation');
    figure(325)
    subplot(1,numel(timescales),ll)
    imagesc(1-squareform(condition_sequence_similarity));
    colorbar;
    title('condition similarity based on sequences')
    %   caxis([-1 1])
    set(gca,'XTick',1:numel(annotation_choose),'XTickLabel',analysisstruct.conditionnames(annotation_choose))
    set(gca,'YTick',1:numel(annotation_choose),'YTickLabel',analysisstruct.conditionnames(annotation_choose))
    xtickangle(90)
    
    %% look at the similarities/diffs across conditions
    if numel(condition_sequence_similarity)
        squared_cond_sim = 1-squareform(condition_sequence_similarity);
        squared_cond_sim_win_across = zeros(2,numel(unique(annotation_categories)));
        squared_cond_sim_win_between = zeros(2,numel(unique(annotation_categories)));
        
        valshere = cell(1,numel(unique(annotation_categories))+1);
        
        for nnk = 1:numel(unique(annotation_categories))
            
            valshere{nnk} = squared_cond_sim(find(annotation_categories==nnk),find(annotation_categories==nnk));
            valshere{nnk} = (triu(valshere{nnk}));
            valshere{nnk}(    valshere{nnk} ==0) = [];
            %  end
            squared_cond_sim_win_across(1,nnk) = nanmean(valshere{nnk});
            squared_cond_sim_win_across(2,nnk) = nanstd(valshere{nnk})./sqrt(numel(valshere{nnk}));
        end
                        hierarchystruct.condition_similarity_within{ll} = valshere;

        for nnk = 1:numel(unique(annotation_categories))
            nnj = 1;
            valshere{nnk} = squared_cond_sim(find(annotation_categories==nnj),find(annotation_categories==nnk));
            valshere{nnk} = valshere{nnk}(:);
            squared_cond_sim_win_between(1,nnk) = nanmean(valshere{nnk});
            squared_cond_sim_win_between(2,nnk) = nanstd(valshere{nnk})./sqrt(numel(valshere{nnk}));
            
        end
        
                hierarchystruct.condition_similarity_tobaseline{ll} = valshere;

        %% bar plots of similarity
        figure(377)
        subplot(3,numel(timescales),ll)
        %      bar(squared_cond_sim_win_across
        %within condition across animals
        bar(squared_cond_sim_win_across(1,:))
        hold on
        errorbar(1:size(squared_cond_sim_win_across,2),squared_cond_sim_win_across(1,:),...
            squared_cond_sim_win_across(2,:),'k','Linestyle','none')
        set(gca,'XTick',1:size(squared_cond_sim_win_across,2),'XTickLabels',...
            nameshere(find(eye(numel(unique(annotation_categories)),numel(unique(annotation_categories))))))
        ylabel('Correlation in sequence structure')
        xtickangle(90)
        box off
       
        
        aggvals = cat(2,squared_cond_sim_win_across(1,:),squared_cond_sim_win_between(1,2:end));
        aggerrs = cat(2,squared_cond_sim_win_across(2,:),squared_cond_sim_win_between(2,2:end));
        aggnames = cat(1,nameshere(find(eye(numel(unique(annotation_categories))))),reshape(nameshere(1,2:end),[],1));
        
        %% dot product between state vectors across each pair of rats/days
        figure(377)
        subplot(3,numel(timescales),numel(timescales)+ll)
        %      bar(squared_cond_sim_win_across
        bar(squared_cond_sim_win_between(1,:),'k')
        hold on
        errorbar(1:size(squared_cond_sim_win_between,2),squared_cond_sim_win_between(1,:),...
            squared_cond_sim_win_between(2,:),'k','Linestyle','none')
        ylabel('Similarity across conditions')
        set(gca,'FontSize',18)
        xtickangle(45)
        set(gca,'XTick',1:size(squared_cond_sim_win_between,2),'XTickLabels',nameshere(1,1:end))
        box off
        
         if ll==numel(timescales)
            print('-depsc',strcat(analysisstruct.plotdirectory,filesep,'Bar_sequence_similarity_USETHIS.epsc'))
            print('-dpng',strcat(analysisstruct.plotdirectory,filesep,'Bar_sequence_similarity_USETHIS.png'))
         end
        
                % hierarchystruct.squared_cond_sim_win_across{ll} = squared_cond_sim_win_across;

        figure(3447)
        subplot(1,numel(timescales),ll)
        bar(aggvals(1,:),'k')
        hold on
        errorbar(1:size(aggvals,2),aggvals,...
            aggerrs,'k','Linestyle','none')
        ylabel('Similarity across conditions')
        set(gca,'FontSize',18)
        xtickangle(45)
        set(gca,'XTick',1:size(aggnames,1),'XTickLabels',aggnames)
        box off
        
        xtickangle(90)
        if ll==numel(timescales)
            print('-depsc',strcat(analysisstruct.plotdirectory,filesep,'Bar_sequence_similarity_var.epsc'))
            print('-dpng',strcat(analysisstruct.plotdirectory,filesep,'Bar_sequence_similarity_var.png'))
        end
        
        % amount of time unstructured
        for kk=setxor(1:numel(annotation_choose),base_condition)
            diff_in_behavior = abs(abs(behavior_vector(kk,:)-behavior_vector(base_condition,:))./(0.5*(behavior_vector(base_condition,:)+behavior_vector(kk,:)))-1);
            [n,x] = ecdf(diff_in_behavior(~isnan(diff_in_behavior)));
            figure(889)
            subplot(1,numel(timescales),ll)
            hold on
            plot(x,n)
            xlabel('dissimilarity in behavioral frequency')
            ylabel('fraction of clusters')
        end
        hold off
        title('cluster frequency rel. baseline')
        legend(analysisstruct.conditionnames(annotation_choose(setxor(1:numel(annotation_choose),base_condition))))
        
    end
    
    
    %% ------------- single rat analyses
    %% get the number of clusters per behavior
    boutlist = [];
    for rr = 1:numel(behs_in_cluster)
        boutlist = cat(1,boutlist,rr*ones(numel(behs_in_cluster{rr}),1));
    end
    allbehshere = unique(base_annotation);
    boutnums = zeros(1,numel(allbehshere));
    for llk = reshape(allbehshere,1,[])
        aahere = cat(2,behs_in_cluster{:});
        boutnums(find(allbehshere==llk)) = numel(unique(boutlist(find((aahere==llk)))));
    end
    
    if numel(boutnums(boutnums>0))
    figure(8778)
    subplot(1,numel(timescales),ll)
    [x,n] = ecdf(boutnums(boutnums>0));
    plot(n,x)
    hold on
    xlabel('number of clusters per behavior')
    ylabel('fraction of behaviors')
    box off
    hold off
    
    
    figure(3993)
    subplot(1,numel(timescales)+2,ll)
    plot_shaded_tsne(analysisstruct,boutnums(2:end),[],1,h1)
    title('number of sequences per behavior')
    hierarchystruct.boutnumsperstate{ll} = boutnums;
    hierarchystruct.behaviorsconsidered = allbehshere;
    
            clustnums =zeros(1,analysisstruct.density_objects);
clustnums(    hierarchystruct.behaviorsconsidered(    hierarchystruct.behaviorsconsidered>0)) = boutnums(hierarchystruct.behaviorsconsidered>0);
h1=figure(5909)
subplot(2,1,1)
        plot_shaded_tsne(analysisstruct,clustnums,[],1,h1)
           title('states most reused')
            
            subplot(2,1,2)
            equivclass_fraction_here =plot_bar_graphs_categorieschanged(analysisstruct,clustnums);
    end
    
    
    %% also look at the overall frequency
    if ll==1
        [un_el,~,ic] = unique(base_annotation);
        counts = accumarray(ic,1);
        beh_shades_here = counts./sum(counts);
        beh_shades_here(1) = [];
        
        subplot(1,numel(timescales)+2,numel(timescales)+1)
        plot_shaded_tsne(analysisstruct,beh_shades_here,[],1,h1)
        title('behavioral frequency')
        
        dwell_times =zeros(1,numel(unique(base_annotation)));
        for mm = reshape(unique(base_annotation),1,[])
            aa = zeros(1,numel(base_annotation));
            aa(find(base_annotation==mm)) = 1;
            rr = bwconncomp(aa);
            dwell_times(mm+1) = mean(cellfun(@numel,rr.PixelIdxList));
        end
        dwell_times(1) = [];
        %dwell_times(dwell_times>15) = 15;
        subplot(1,numel(timescales)+2,numel(timescales)+2)
        plot_shaded_tsne(analysisstruct,dwell_times,[],1,h1)
        title('dwell times')
    end
    
    [nn,xx] = ecdf(cellfun(@numel,conditions_using_seq));
    figure(8989)
    subplot(1,numel(timescales),ll)
    plot(xx,nn)
    xlabel('number of conditions')
    ylabel('fraction of sequences')
    
    figure(9034)
    subplot(1,numel(timescales),ll)
    imagesc(1-squareform(pdist(base_annot_vectors,'correlation')))
    set(gca,'XTick',1:numel(annotation_choose),'XTickLabels',analysisstruct.conditionnames(annotation_choose));
    set(gca,'YTick',1:numel(annotation_choose),'YTickLabels',analysisstruct.conditionnames(annotation_choose));
    xtickangle(90)
    colorbar;
    
    %% also get the distribution of
    %[nnmf_annotation{ll},nnmf_score{ll},D] =nnmf(annotation_matrix_shifted_conv_agg{ll},numnnmf_pcs(ll));
    %
    %     if strcmp(cases{kkjj},'real')
    %         nnmf_score_out{ll} = nnmf_score{ll};
    %     end
    %
    %     annot_resampled= annotation_matrix_shifted_conv_agg{ll}(:,1:decimation_factor:end);
    %     %Z = linkage(annot_resampled','average','correlation');
    %     Z = linkage(annot_resampled','ward','euclidean');
    %
    %     %T = cluster(Z,'cutoff',)
    % %     figure(4)
    %     if sort_dendrogram
    %         [H,T,outperm] =dendrogram(Z,size(annotation_matrix_shifted_conv_agg{ll}(:,1:decimation_factor:end),2));
    %         %xtickangle(90)
    %         outperms{ll} = outperm;
    %     end
    % c = cluster(Z,'maxclust',15);
    % [newind,sorted_clust_ind] = sort(c,'ASCEND');
    
end

hierarchystruct.clustered_behavior = clustered_behavior;
hierarchystruct.equivclass_fraction_agg = equivclass_fraction_agg;
hierarchystruct.behavior_composition_vector = behavior_composition_vector;
hierarchystruct.normshades_agg = normshades_agg;
hierarchystruct.normshades_agg = normshades_agg;
hierarchystruct.ratcond_names = ratcond_names;%analysisstruct.conditionnames(ratcond_inds);
hierarchystruct.annotation_choose = annotation_choose;%([1,2,3,7]);
hierarchystruct.base_annotation = base_annotation;
hierarchystruct.base_annotation_inds = base_annotation_inds;
hierarchystruct.annotation_indexes=annotation_indexes;
hierarchystruct.timescales=timescales;
hierarchystruct.uniquevals=uniquevals;
hierarchystruct.timescales=timescales;

%Params, for ASD, caff, etc
%return
%                   make_layer_cake_from_hierarchicalannotation(hierarchystruct.clustered_behavior);

%make_state_example_plots(analysisstruct,hierarchystruct,timescales,uniquevals)




%% only do detailed if not looking at coarse for now
if params.do_cluster_analysis
    if ~strcmp(cases{kkjj},'coarse')
        
        %% Mutual Information between seq/state and state/behavior
        MutualInformation(clustered_behavior{2},clustered_behavior{1})
        MutualInformation(clustered_behavior{2},base_annotation)
        MutualInformation(clustered_behavior{2},clustered_behavior{1})
        annotation_choose
        
        %% look at the influence of behavioral state on this timescale
        circshiftval =10; %non decimated, for beh similarity
        cluster_compare = 2; % for behavioral similarity
        
        
        
        %
        %% explore tsne
        h1=figure(607)
        params.nameplot=1;
        params.density_plot =0;
        params.watershed = 1;
        params.sorted = 1;
        params.do_coarse = 0;
        
        plot_clustercolored_tsne(analysisstruct,1,params.watershed,h1,params)
        set(gcf,'renderer','opengl')
        axis off
        
        %% get the annotation matrix on some timescale
        annotation_use = base_annotation;
        
        annotation_matrix_shifted = zeros(numel(unique(base_annotation)),numel(base_annotation));
        for rr = 1:(max(base_annotation))
            annotation_matrix_shifted(rr,find(annotation_use==rr)) = 1;
        end
        %annotation_matrix_shifted(find(sum(annotation_matrix_shifted,2)==0),:) =[];
        %else
        %    annotation_matrix_shifted = nnmf_score{ll-1};
        %end
        
        % get PDF of behaviors at each timepoint for chosen timescale
        annotation_matrix_shifted_conv = conv2(annotation_matrix_shifted,ones(1,circshiftval),'same');
        annotation_matrix_shifted_conv_agg = bsxfun(@rdivide,annotation_matrix_shifted_conv,nansum(annotation_matrix_shifted_conv,1));
        annotation_matrix_shifted_conv_agg(isnan(annotation_matrix_shifted_conv_agg)) = 0;
        
        %for each behavioral state get the behavioral similarity -- how similar
        %is the behavioral mileau on different clustering instances. How
        %stereotyped are different epochs
        beh_similarity = cell(1,numel(unique(clustered_behavior{cluster_compare})));
        for kk= unique(clustered_behavior{cluster_compare})
            beh_similarity{1+kk} = 1-squareform(pdist(annotation_matrix_shifted_conv_agg(:,find(clustered_behavior{cluster_compare}==kk)),'correlation'));
            figure(99)
            subplot(1,numel(unique(clustered_behavior{cluster_compare})),1+kk)
            imagesc(beh_similarity{kk+1})
        end
        
        
        %% assemble the transitions for each behavior
        numshuff = 10;
        per_behavior_similarities = cell(1,analysisstruct.density_objects);
        similarity_difference = zeros(1,analysisstruct.density_objects);
        similarity_difference_full = cell(1,analysisstruct.density_objects);
        similarity_difference_mostmodpair = cell(1,analysisstruct.density_objects);
        condentropy_sum =zeros(2,analysisstruct.density_objects);
        featureval_states_beh = nan*zeros(analysisstruct.density_objects,numel(unique(clustered_behavior{cluster_compare})),size(analysisstruct.jt_features,2));
        featureval_states_beh_seq = nan*zeros(analysisstruct.density_objects,numel(unique(clustered_behavior{cluster_compare})),...
            numel(unique(clustered_behavior{1})),size(analysisstruct.jt_features,2));
         featureval_states_beh_seq_random = nan*zeros(analysisstruct.density_objects,numel(unique(clustered_behavior{cluster_compare})),...
            numel(unique(clustered_behavior{1})),size(analysisstruct.jt_features,2),numshuff);
        
        instances_states_beh = zeros(analysisstruct.density_objects,numel(unique(clustered_behavior{cluster_compare})));
       % normalize the features to enable comparison
        normjtfeat = bsxfun(@rdivide,bsxfun(@minus,analysisstruct.jt_features,nanmean(analysisstruct.jt_features,1)),nanstd(analysisstruct.jt_features,[],1));
        normjtfeat(normjtfeat>5) = 5;
        normjtfeat(normjtfeat<-5) = -5;
        featureval_av = zeros(analysisstruct.density_objects,size(analysisstruct.jt_features,2));
              featureval_std = zeros(analysisstruct.density_objects,size(analysisstruct.jt_features,2));
  
        %number of instances to keep for analyses
        inst_thresh = 40;
        seq_inst_thresh = 5;
        
        do_cond_entropy = 0;
        if (do_cond_entropy)
            for nn = setxor(unique(base_annotation),intersect(unique(base_annotation),0))
                %% get the total conditional entropy for each behavior in context explaining its sequential mileau H(seq|state)./H(seq). A perfect value would mean that
                %there were two sequences totally explained by the state for ea
                %behavior
                %problem is N2 calls
                for jj = setxor(unique(base_annotation),intersect(unique(base_annotation),0))
                    if numel(find(base_annotation==nn))
                        if condentropy(  annotation_matrix_shifted_conv_agg(jj,find(base_annotation==nn)))
                            condentropy_sum(1,nn)= condentropy_sum(1,nn)+ condentropy(  annotation_matrix_shifted_conv_agg(jj,find(base_annotation==nn)),...
                                clustered_behavior{cluster_compare}(find(base_annotation==nn)))...
                                ./condentropy(  annotation_matrix_shifted_conv_agg(jj,find(base_annotation==nn)));  %how variable is the sequencing
                        end
                    end
                    %condentropy_sum(1,nn)= condentropy_sum(1,nn)+condentropy(  annotation_matrix_shifted_conv_agg(jj,find(base_annotation==nn)));
                end
            end
        end
        
        
        
        
        for nn = reshape(setxor(unique(base_annotation),intersect(unique(base_annotation),0)),1,[])
            %% look at the similarities -- for each behavior how similar is its mileau on one timescale given the state label
            per_behavior_similarities{nn} = zeros( numel(unique(clustered_behavior{cluster_compare})),max(base_annotation));
            for kk= unique(clustered_behavior{cluster_compare})
                per_behavior_similarities{nn}(kk+1,2:end)  = beh_similarity{kk+1}(nn,2:end);
                
                if numel(intersect(find(clustered_behavior{cluster_compare} == kk),find(base_annotation==nn))) <inst_thresh
                    per_behavior_similarities{nn}(kk+1,2:end) = zeros(1,numel(beh_similarity{kk+1}(nn,2:end)));
                end
            end
        end
        
        
        
        per_behavior_similarities{nn}(:,nn) = 0;
        per_behavior_similarities{nn}(isnan(per_behavior_similarities{nn})) = 0;
        similarity_difference(nn) = nanmean(1-pdist(per_behavior_similarities{nn},'correlation'));
        similarity_difference_full{nn} = squareform(1-pdist(per_behavior_similarities{nn},'correlation'));
        minval = min(similarity_difference_full{nn}(similarity_difference_full{nn}>0));
        if numel(minval)
            [ind1 ind2] = ind2sub(size(similarity_difference_full{nn}),find(similarity_difference_full{nn}==minval,1,'first'));
            similarity_difference_mostmodpair{nn} = [nn minval ind1 ind2];
        else
            similarity_difference_mostmodpair{nn} = [0 0 0 0 ];
        end
        
        
        similarity_difference_mostmodpair_mat= cell2mat(similarity_difference_mostmodpair');
        
        
        [sortval,sortind] = sort(similarity_difference_mostmodpair_mat(:,2),'DESCEND');
        similarity_difference_mostmodpair_mat(sortind(sortval>0),:)
        
        
        
        %% get the average jt features, to see which are most modulated
                for nn = reshape(setxor(unique(base_annotation),intersect(unique(base_annotation),0)),1,[])
                              behavior_ind = reshape(find(base_annotation==nn),[],1);
                              behavior_ind(behavior_ind>size(normjtfeat,1)) = size(normjtfeat,1);
            for kk= unique(clustered_behavior{cluster_compare})
                        instances_states_beh(nn,kk+1) = numel(behavior_ind);
            end
            featureval_av(nn,:) = squeeze(nanmean(normjtfeat(behavior_ind,:),1));
            featureval_std(nn,:) = squeeze(nanstd(normjtfeat(behavior_ind,:),[],1));
                end
        
           
      
              
              %% hash these
              loopseq = setxor(unique(clustered_behavior{1}),0);
              for mmk = loopseq
                  sequence_ind_agg{mmk} = find(clustered_behavior{1} == mmk);
              end
              loopstate = setxor(unique(clustered_behavior{cluster_compare}),0);
              
              for kk= loopstate
                  state_ind_agg{kk} = find(clustered_behavior{cluster_compare} == kk);
              end
              
              
              %% loop over base annotation
              tic
              for nn = setxor(unique(base_annotation),...
                      intersect(unique(base_annotation),0))                 
                  beh_here = find(base_annotation==nn);
                  if numel(beh_here>inst_thresh)
                      for mmk = loopseq
                          seq_beh_overlap = intersect(beh_here,sequence_ind_agg{mmk});
                          if numel(seq_beh_overlap)>seq_inst_thresh
                              for kk= loopstate
                                  behavior_ind = intersect(seq_beh_overlap,state_ind_agg{kk});
                                  if numel(behavior_ind) >seq_inst_thresh
                                      if mmk == loopseq(1)
                                      featureval_states_beh(nn,kk+1,:) = nanmean(normjtfeat(intersect(state_ind_agg{kk},beh_here),:),1);
                                      end
                                      featureval_states_beh_seq(nn,kk+1,mmk+1,:) = nanmean(normjtfeat(behavior_ind,:),1);
                                      for ll=1:numshuff
                                          featureval_states_beh_seq_random(nn,kk+1,mmk+1,:,ll) = nanmean(normjtfeat(...
                                              randsample(seq_beh_overlap,numel(behavior_ind)),:),1);                                      
                                  end
                              end
                          end
                      end
                  end
                  end
              end
              toc
            
                          
              
        for nn = setxor(unique(base_annotation),intersect(unique(base_annotation),0))            
            for kk= unique(clustered_behavior{cluster_compare})
                behavior_ind = intersect(find(clustered_behavior{cluster_compare} == kk),find(base_annotation==nn));
                behavior_ind(behavior_ind>size(normjtfeat,1)) = size(normjtfeat,1);
              %  if numel(behavior_ind) >inst_thresh
                    featureval_states_beh(nn,kk+1,:) = nanmean(normjtfeat(behavior_ind,:),1);
                   % for mmk = unique(clustered_behavior{1})
                   %     behavior_ind_seq = intersect(find(clustered_behavior{1} == mmk),behavior_ind);
                   %     if numel(behavior_ind_seq>seq_inst_thresh)
                   %         featureval_states_beh_seq(nn,kk+1,mmk+1,:) = nanmean(normjtfeat(behavior_ind_seq,:),1);
                   %     end
                   % end
                end
            end
        
%         %% find a juicy cluster to investigate
%         h2 = figure(89)
%         plot(analysisstruct.zValues(:,1),analysisstruct.zValues(:,2),'+')
%         buffer = [-1:1]
%         % beh_list = examine_features(h2,analysisstruct.zValues,...
%         %     analysisstruct.subset_of_points_to_plot_tsne_capped,...
%         %     analysisstruct.condition_inds,  analysisstruct.subset_of_points_to_plot_tsne_capped,...
%         %     analysisstruct.mocapstruct_reduced_agg,[0],1:numel(analysisstruct.conditionnames),analysisstruct.conditionnames,buffer);
%         clusterpicked = 12;% mode(analysisstruct.annot_reordered{1}(beh_list{1}));
%         indpick = 1;%find(sortind == mode(analysisstruct.annot_reordered{1}(beh_list{1})));
        %animate_markers_aligned_fullmovie(analysisstruct.mocapstruct_reduced_agg{1},intersect(find(clustered_behavior{3}==5),find(base_annotation==55)))
        
        %% compute modulation over states controlling for sequences
        %zscoreddiffs_seqcont = squeeze(nanstd(featureval_states_beh_seq,[],2)./nanmean(featureval_states_beh_seq,2));
        shuff_random = squeeze(nanmean(featureval_states_beh_seq_random,5));
        fprintf('looking at jt feature change controlling for sequence -- aggregating\n')
        diff_overstate_seqcont = nan*zeros(size(featureval_states_beh_seq,2),size(featureval_states_beh_seq,2),size(featureval_states_beh_seq,1) );
        diff_overstate_seqcont_random = ...
            nan*zeros(size(featureval_states_beh_seq,2),size(featureval_states_beh_seq,2),size(featureval_states_beh_seq,1),numshuff );
tic
featdiff_agg = [];
featdiff_agg_shuff = [];

        for jk = 1:size(featureval_states_beh_seq,2)
            for jkj = 1:size(featureval_states_beh_seq,2)
                %look at states with the same sequences
                goodseq_1= find(nansum(nansum(featureval_states_beh_seq(:,jk,:,:),1),4));
                goodseq_2= find(nansum(nansum(featureval_states_beh_seq(:,jkj,:,:),1),4));
                goodseq = intersect(goodseq_1,goodseq_2);
                %compute the difference (note norm jt features is already
                %z-scored)
                feat_diff = squeeze(nanmean(featureval_states_beh_seq(:,jkj,goodseq,:)-featureval_states_beh_seq(:,jk,goodseq,:),3));
               
                valstosave = reshape(nanmean(bsxfun(@rdivide,...
                    abs(squeeze(featureval_states_beh_seq(:,jkj,goodseq,:)-...
                    featureval_states_beh_seq(:,jk,goodseq,:))),reshape(featureval_std,size(featureval_std,1),1,size(featureval_std,2))),2),[],1);
             valstosave(valstosave==0) = [];
                featdiff_agg = cat(1,featdiff_agg,valstosave(find(~isnan(valstosave))));
                %fractional change in mean in the units of feature standard
                %deviations
                feat_diff_norm = abs(feat_diff)./(featureval_std);
                diff_overstate_seqcont(jk,jkj,:) = nanmean(feat_diff_norm,2);
                
                %also shuffle the state labels 
                for ll=1:numshuff
            feat_diff_random = squeeze(nanmean(featureval_states_beh_seq_random(:,jkj,goodseq,:,ll)-featureval_states_beh_seq_random(:,jk,goodseq,:,ll),3));

                diff_overstate_seqcont_random(jk,jkj,:,ll) = nanmean(abs(feat_diff_random)./(featureval_std),2);
                
                valstosave = reshape(nanmean(bsxfun(@rdivide,...
                    abs(squeeze(featureval_states_beh_seq_random(:,jkj,goodseq,:,ll)-...
                    featureval_states_beh_seq_random(:,jk,goodseq,:,ll))),...
                    reshape(featureval_std,size(featureval_std,1),1,size(featureval_std,2))),2),[],1);
            
                valstosave(valstosave==0) = [];
                                featdiff_agg_shuff = cat(1,featdiff_agg_shuff,valstosave(find(~isnan(valstosave))));

                
                end
            end
        end
        toc
        
        %get nonzero state comparisons
        diff_overstate_seqcont(diff_overstate_seqcont==0) = nan; %remove diagonal
                diff_overstate_seqcont_random(diff_overstate_seqcont==0) = nan; %remove diagonal

        % look at average feature change over all pairs of states
        val_diffs = squeeze(nanmean(nanmean(diff_overstate_seqcont,1),2));
        val_diffs(isnan(val_diffs)) = 0;
        
         val_diffs_random = squeeze(nanmean(nanmean(nanmean(diff_overstate_seqcont_random,1),4),2));
        val_diffs_random(isnan(val_diffs_random)) = 0;
        
        val_diffs(find(val_diffs))
        val_diffs_random(find(val_diffs_random))
        
        % look at aggregate of all behavior-sequence pairs
        featdiff_agg_shuff(isnan(featdiff_agg_shuff)) = 0;
        
        
        bins = linspace(0,prctile(featdiff_agg,99),50);
        figure(8887)
       %[yvals,n] = hist(val_diffs(val_diffs>0),bins);
              [yvals,n] = hist(featdiff_agg(featdiff_agg>0),bins);
       plot(bins,yvals./sum(yvals),'r');
        hold on
                      [yvals_n,n] = hist(featdiff_agg_shuff(featdiff_agg_shuff>0),bins);
          %   [yvals_r,n] = hist(val_diffs_random(val_diffs_random>0),bins);
         plot(bins,yvals_n./sum(yvals_n),'k');

hold off
   hierarchystruct.featdiff_agg_shuff = featdiff_agg_shuff;
   hierarchystruct.featdiff_agg = featdiff_agg;


        hh=figure(8888)
        subplot(2,1,1)
        plot_shaded_tsne_absolute_alpha(analysisstruct,val_diffs,[],1,h1)
        title('controlling for sequential context, states with most state-kin coupling (av diff/mean)')
        axis off
        hierarchystruct.statemodulation_of_kinematics = val_diffs;
            hierarchystruct.statemodulation_of_kinematics_shuffled = val_diffs_random;
    
        val_diffs(val_diffs==0) = nan;
        subplot(2,1,2)
        equivclass_fraction_here =plot_bar_graphs_categorieschanged(analysisstruct,val_diffs);
        if isfield(analysisstruct,'behavior_order')
            bar(equivclass_fraction_here(analysisstruct.behavior_order_numbers) ,'FaceColor',[0.5 .5 .5],'EdgeColor','none')
            set(gca,'XTickLabel',analysisstruct.behavior_order)
            xtickangle(90)
        end
        hh.Renderer='Painters';
        print('-depsc',strcat(analysisstruct.plotdirectory,filesep,'state_kin_couple.epsc'))
        print('-dpng',strcat(analysisstruct.plotdirectory,filesep,'state_kin_couple.png'))
        
        hierarchystruct.state_kin_coupling_contseq = val_diffs;
        
        %
        % zscoreddiffs_seqcont = featureval_states_beh_seq(:,jk,:,:)
        % zscoreddiffs = abs(zscoreddiffs);
        % zscoreddiffs(zscoreddiffs==0) = nan;
        % zscoreddiffs(zscoreddiffs>5) = 5;
        
        
        %% compute feature modulation so see which are most modulated across states
        zscoreddiffs = squeeze(nanstd(featureval_states_beh,[],2)./nanmean(featureval_states_beh,2));
        zscoreddiffs = abs(zscoreddiffs);
        zscoreddiffs(zscoreddiffs==0) = nan;
        zscoreddiffs(zscoreddiffs>5) = 5;
        modulation_index = squeeze(nanmean(zscoreddiffs,2));
        modulation_index(isnan(modulation_index)) = 0;
        
        %also normalize the
        featureval_av = abs(featureval_av);featureval_av(featureval_av>5) = 5;
        featureval_av = squeeze(nanmean(featureval_av,2));
        
        hh=figure(33)
        subplot(1,3,1)
        plot(modulation_index)
        title('modulation index of jt features (av z in beh)')
        %modulation_index(modulation_index>2) =2;
        subplot(1,3,2)
        plot_shaded_tsne_absolute(analysisstruct,modulation_index./featureval_av,[],1,h1)
        title('modulation index per beh')
        
        %plot_shaded_tsne(analysisstruct,1:152,[],1,h1)
        subplot(1,3,3)
        plot_shaded_tsne(analysisstruct,featureval_av,[],1,h1)
        title('av feature value')
        hh.Renderer='Painters';
        print('-depsc',strcat(analysisstruct.plotdirectory,filesep,'overallmodulationacrossstates.epsc'))
        print('-dpng',strcat(analysisstruct.plotdirectory,filesep,'overallmodulationacrossstates.png'))
        
        hierarchystruct.statemodulation_of_kinematics_altmetric = modulation_index./featureval_av;
        
        
       
        
        
        
        %featuremat(find(isnan(sum(featuremat,2))),:) = [];
        %   framesuse = unique(  bsxfun(@plus, intersect(find(clustered_behavior{cluster_compare} == badind1-1),find(base_annotation==clusterpicked)),(-3:1:3)'));
        %animate_markers_aligned_fullmovie(analysisstruct.mocapstruct_reduced_agg{1},framesuse)
        %   framesuse = unique(  bsxfun(@plus, intersect(find(clustered_behavior{cluster_compare} == badind2-1),find(base_annotation==clusterpicked)),(-3:1:3)'));
        %animate_markers_aligned_fullmovie(analysisstruct.mocapstruct_reduced_agg{1},framesuse)
        
        plotfolder = strcat('Y:\Jesse\Data\Motionanalysis_captures\Cluster_videos\20190205_state_behavior_sequencing\');
        mkdir(plotfolder)
        dovideo = 0;
        indrel = cell(1,2);
        clusterpicked=13;
         figure(35)
        %subplot(1,3,3)
        imagesc(squeeze(featureval_states_beh(clusterpicked,:,:)))
        title('feature values across states')
        %% Make example plots and videos comparing example state with biggest discrepancy
        featuremat = squeeze(featureval_states_beh(clusterpicked,:,:));
        statediff = 1-squareform(pdist(featuremat,'correlation'));
        [valmin] = min(statediff(:));
        [badind{1} badind{2}] = ind2sub(size(featuremat),find(statediff==valmin,1,'first'));
        
        if isfield(analysisstruct,'agg_mocap_structs_snippets') && numel(analysisstruct.agg_mocap_structs_snippets)
            for mii =1:2
                behvals =  intersect(find(clustered_behavior{cluster_compare} == badind{mii}-1),find(base_annotation==clusterpicked));
                indrel{mii} = find(ismember(analysisstruct.agg_snippetinds{1},behvals(1:5:end) ));
                mocapstruct_fullres = analysisstruct.mocapstruct_reduced_agg{1};
                mocapstruct_fullres.markers_aligned_preproc = analysisstruct.agg_mocap_structs_snippets{1};
                mocapstruct_fullres.markers_preproc = analysisstruct.agg_preproc{1};
                % mocapstruct_fullres.markers_aligned_preproc = agg_mocap_preproc;
                if dovideo
                    
                    animate_markers_aligned_fullmovie(mocapstruct_fullres,indrel{mii}(1:4:end))
                    
                    v = VideoWriter(strcat(plotfolder,filesep,'clusteredmovie_',num2str(clusterchoose),'_state_',num2str((statechoose))),'MPEG-4');
                    open(v)
                    writeVideo(v, M)
                    close(v)
                end
                % animate_markers_aligned_fullmovie(analysisstruct.mocapstruct_reduced_agg{1},behvals)
                
                %   indrel2 = find(ismember(analysisstruct.agg_snippetinds{1},...
                %   intersect(find(clustered_behavior{cluster_compare} == badind2-1),find(base_annotation==clusterpicked))));
                figure(101)
                subplot(1,2,mii)
                
                plot(analysisstruct.agg_mocap_structs_snippets{1}.HeadF(indrel{mii},:))
                box off
                xlabel('frame (at ~60Hz)')
                ylabel('Trunk Position')
            end
            mocapstruct_fullres.plotdirectory = [];
            %  compare_plot_marker_characteristics_timerange(mocapstruct_fullres,indrel{1},mocapstruct_fullres,indrel{2})
        end
        
        %% plot the dissimilarities in sequence ID
        figure(904)
        subplot(1,2,1)
        bar(sort(1-similarity_difference,'DESCEND'))
        title('average dissimilarity in sequencing across all clusters')
        subplot(1,2,2)
        bar(1-similarity_difference_mostmodpair_mat(sortind,2))
        title('peak dissimilarity in sequencing')
        
        hierarchystruct.analog_seqvariability = similarity_difference;
        
        
        figure(304)
           subplot(2,1,1)
            plot_shaded_tsne_absolute(analysisstruct,condentropy_sum(1,:),[],1,h1)
            title('states most reused')
            
            subplot(2,1,2)
            equivclass_fraction_here =plot_bar_graphs_categorieschanged(analysisstruct,condentropy_sum(1,:));
            
        
        % axis is the normalized (ratio) conditional entropy
        if (do_cond_entropy)
            
            figure(303)
            subplot(2,1,1)
            plot_shaded_tsne_absolute(analysisstruct,condentropy_sum(1,:),[],1,h1)
            title('states most reused')
            
            subplot(2,1,2)
            equivclass_fraction_here =plot_bar_graphs_categorieschanged(analysisstruct,condentropy_sum(1,:));
            
            hh.Renderer='Painters';
            print('-depsc',strcat(analysisstruct.plotdirectory,filesep,'state_behreuse',num2str(cluster_compare),'.epsc'))
            print('-dpng',strcat(analysisstruct.plotdirectory,filesep,'state_behreuse',num2str(cluster_compare),'png'))
            
        end
        %  condentropy_sum
        
        
        %% get the cluster ID of the cluster of choice
        
        %% pick the cluster
        indpick= 0;
        if indpick>0
        clusterchoose = similarity_difference_mostmodpair_mat(sortind(indpick),1);
        states_here = [badind{1} badind{2}];
        %states_here = similarity_difference_mostmodpair_mat(sortind(indpick),3:4);
        figure(809)
        %imagesc(per_behavior_similarities{similarity_difference_mostmodpair_mat(sortind(indpick),1)}(states_here ,:))
        %figure(809)
        %imagesc(per_behavior_similarities{similarity_difference_mostmodpair_mat(sortind(indpick),1)})
        %% watch a behavior in two different contexts
        
        % also 54 and 3 and 5
        dovideo=0;
        for statechoose = states_here
            %statechoose = 3;
            framecluster = intersect(find(clustered_behavior{cluster_compare} == statechoose-1),find(base_annotation==clusterchoose));
            beh_shades = zeros(1,numel(unique(base_annotation)));
            [un_el,~,ic] = unique(base_annotation(framecluster));
            counts = accumarray(ic,1);
            beh_shades(un_el+1) = counts./sum(counts);
            beh_shades(1) = [];
            
            figure(301)
            subplot(1,numel(states_here)+1,1)
            plot_shaded_tsne(analysisstruct,beh_shades,[],1,h1)
            
            %% get the values for the sequential context
            framesuse = unique(  bsxfun(@plus, intersect(find(clustered_behavior{cluster_compare} == statechoose-1),...
                find(base_annotation==clusterchoose)),(-20:1:20)'));
            framesuse(framesuse<1) = 1;
            framesuse(framesuse>numel(base_annotation)) = numel(base_annotation);
            
            beh_shades = zeros(1,numel(unique(base_annotation)));
            [un_el,~,ic] = unique(base_annotation(framesuse));
            counts = accumarray(ic,1);
            beh_shades(un_el+1) = counts./sum(counts);
            beh_shades(1) = [];
            
            figure(301)
            subplot(1,numel(states_here)+1,find(statechoose==states_here)+1)
            plot_shaded_tsne(analysisstruct,beh_shades,[],1,h1)
            
            if dovideo
                M=animate_markers_nonaligned_fullmovie(analysisstruct.mocapstruct_reduced_agg{1},framesuse(1:min(400,numel(framesuse))))
                v = VideoWriter(strcat(plotfolder,filesep,'clusteredmovie_',num2str(clusterchoose),'_state_',num2str((statechoose))),'MPEG-4');
                open(v)
                writeVideo(v, M)
                close(v)
            end
        end
        if dovideo
            framesuse = unique(  bsxfun(@plus, find(base_annotation==clusterchoose),(-1:1:1)'));
            M=animate_markers_nonaligned_fullmovie(analysisstruct.mocapstruct_reduced_agg{1},framesuse(1:min(400,numel(framesuse))))
            v = VideoWriter(strcat(plotfolder,filesep,'clusteredmovie_',num2str(clusterchoose)),'MPEG-4');
            open(v)
            writeVideo(v, M)
            close(v)
        end
        
        
        % M=animate_markers_nonaligned_fullmovie_white(analysisstruct.mocapstruct_reduced_agg{1},20000:3:30000)
        
        %% look at the same b
        beh_similarity_all = zeros(1,numel(unique(base_annotation)));
        
        for kk= unique(base_annotation)
            indshere = find(base_annotation==kk);
            beh_similarity_all(1+kk) = nanmean(1-(pdist(annotation_matrix_shifted_conv_agg(:,indshere(1:5:end))','correlation')));
        end
        beh_similarity_all(1) = [];
        %% get shaded behavior
        beh_shades = beh_similarity_all;
        
        examine_feat_here = 0;
        if examine_feat_here
            figure(302)
            plot_shaded_tsne(analysisstruct,beh_shades,[],1,h1)
            h2 = figure(89)
            plot(analysisstruct.zValues(:,1),analysisstruct.zValues(:,2),'+')
            buffer = [-1:1]
            beh_list = examine_features(h2,analysisstruct.zValues,...
                analysisstruct.subset_of_points_to_plot_tsne_capped,...
                analysisstruct.condition_inds,  analysisstruct.subset_of_points_to_plot_tsne_capped,...
                analysisstruct.mocapstruct_reduced_agg,[1],1:numel(analysisstruct.conditionnames),analysisstruct.conditionnames,buffer);
            
            %% find behaviors with most interaction between kinematics and variable state
            analysisstruct.tsnefeat_name_dotplot
        end
        
        %    figure(199)
        % %subplot(1,numel(unique(clustered_behavior{cluster_compare})),1+kk)
        % imagesc(beh_similarity_all{1+kk})
        
        end
    end
end

end

