%% construct a behavioral hierarchy from an input analyssstruct
function [hierarchystruct]= find_sequences_states_demo(analysisstruct,condname_in,params)

if nargin<3
    params.do_show_pdistmatrix =0;
    params.cases ='real';
    params.decimation_factor = 10;
end

if ~isfield(analysisstruct,'plotdirectory')
    analysisstruct.plotdirectory = '';
end

if ~isfield(params,'timescales')
    timescales= [1./4  2 ];
else
    timescales = params.timescales;
end



if isfield(params,'corr_threshold')
    corr_threshold = params.corr_threshold;
    clustercutoff = params.clustercutoff;
else
    corr_threshold = 0.3;
    clustercutoff = 0.65;
end
chopsize = 20000;



if ~isfield(params,'doclustering')
    params.doclustering = 1;
end


condnames = condname_in;

ratcond_inds = [];
annotation_categories = [];

%% -----------------------------------------
annotation_choose=1;

decimation_factor = params.decimation_factor;
seq_inds_cluster_agg= cell(1,numel(timescales));

%% loop over all elements
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

%% matricies for saving
clustered_behavior = cell(1,numel(timescales));
clustered_sequences = cell(1,numel(timescales));
sequences_inds = cell(1,numel(timescales));

hierarchystruct = struct();
hierarchystruct.annotation_categories = annotation_categories;


for ll = 1:numel(timescales)
    circshiftval = round((300./(analysisstruct.tsnegranularity))*60*timescales(ll));
    annotation_use = base_annotation;
    
    annotation_matrix_shifted = zeros(max((base_annotation)),numel(base_annotation));
    for rr = reshape((unique(base_annotation(base_annotation>0))),1,[])
        annotation_matrix_shifted(rr,find(annotation_use==rr)) = 1;
    end
    annotation_matrix_shifted(find(sum(annotation_matrix_shifted,2)==0),:) =[];
    
    % get PDF of behaviors at each timepoint
    annotation_matrix_shifted_conv = conv2(annotation_matrix_shifted,ones(1,circshiftval),'same');
    annotation_matrix_shifted_conv_agg = bsxfun(@rdivide,annotation_matrix_shifted_conv,nansum(annotation_matrix_shifted_conv,1));
    annotation_matrix_shifted_conv_agg(isnan(annotation_matrix_shifted_conv_agg)) = 0;
    
    
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
                [sub1,sub2] = ind2sub(size(L),density_cc_temp.PixelIdxList{rr});
                density_cc.PixelIdxList{end+1} = sub2ind([size(annotation_matrix_shifted_conv_agg,2) size(annotation_matrix_shifted_conv_agg,2)],...
                    bsxfun(@plus,sub1,(jj-1)*chopsize),bsxfun(@plus,sub2,(kk-1)*chopsize));
                goodobj = goodobj+1;
            end
            density_cc.NumObjects = density_cc.NumObjects+goodobj;
        end
    end
    clear vals_examine_here pdistmatrixhere density_cc_temp
    
    figure(5454)
    subplot(1,numel(timescales),ll)
    bar(corr_hist_range,histsum./sum(histsum))
    box off
    ylim([0 0.2])
    
    
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


fprintf('now starting clustering \n')
equivclass_fraction_agg = cell(1,numel(timescales));
normshades_agg = cell(1,numel(timescales));
behshades_agg = cell(1,numel(timescales));
behavior_composition_vector = cell(1,numel(timescales));

for ll = 1:numel(timescales)
    %% clustering step
    fprintf('correlation based clustering over %f points \n',size(clustered_sequences{ll},2))
    
    
        behsuse = 1:size(clustered_sequences{ll},1);
    
  
        clustercutoff_here = clustercutoff;
    
    
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
    
    %% setup the loop over clusters
    seq_inds_cluster_agg{ll} = cell(1,numel(unique(T)));
    conditions_using_seq =cell(1,numel(unique(T)));
    base_annot_clustered = zeros(1,numel(base_annotation));
    base_annot_vectors = zeros(numel(unique(annotation_choose)),numel(unique(T)));
    
    %% get the number of 'good' clusters
    [un_el,~,ic] = unique(T);
    counts_seq = accumarray(ic,1);
    goodclust = find(counts_seq>3);
    fprintf('numel clust (min 20 fr) %f for timescale %f \n',  numel(goodclust),timescales(ll)*60);
    [~,gcsorted] = sort(counts_seq(goodclust),'DESCEND');
       
    allbehshere = reshape(unique(base_annotation),1,[]);
    do_save_videos=0;
    beh_shades_agg = zeros(numel(gcsorted),analysisstruct.density_objects+1);
    behs_in_cluster = [];
    for clusterhere = 1:numel(gcsorted)
        clusterinds = find(T==goodclust(gcsorted(clusterhere)));

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

        normshades = log2(beh_shades);
        normshades(isinf(normshades)) = -13;
        normshades_agg{ll}{clusterhere} = normshades;
        behshades_agg{clusterhere} = beh_shades_agg;
        

        figure(395+ll)
        h1=subplot_tight(7,3,clusterhere)

        plot_shaded_tsne_absolute_alpha(analysisstruct,normshades,[],1,h1,1)
        axis equal
        axis off
        colorbar off
        caxis([0.01 0.25])
        if strcmp(cases{kkjj},'coarse')
            fprintf('for cluster %f frac move: %f frac fast %f frac rest %f \n',clusterhere, beh_shades(end-2),beh_shades(end-1),beh_shades(end));
        end
        if clusterhere == numel(gcsorted)
            mkdir(analysisstruct.plotdirectory)
            set(gcf,'renderer','opengl')           
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




end

