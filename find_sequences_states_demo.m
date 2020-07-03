function [hierarchystruct]= find_sequences_states_demo(analysisstruct,annotation_choose,params)
%% construct a behavioral hierarchy from an input analyssstruct
% ---------------------------
% (C) Jesse D Marshall 2020
%     Harvard University 




%% load default parameters
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

annotation_categories = [];

%% -----------------------------------------
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
   
    %% break up watershed into blocks to get over memory complexity issue
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


%% cluster the data, and plot the sequences
fprintf('now starting clustering \n')
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
    end

     clustered_behavior{ll} = base_annot_clustered;     

     behavior_vector = zeros(numel(annotation_choose),max(clustered_behavior{ll})+1);
    for kk =1:numel(annotation_choose)
        beh_use_vec =  clustered_behavior{ll}(find(base_annotation_inds ==annotation_choose(kk)));
        [un_el,~,ic] = unique(beh_use_vec);
        counts = accumarray(ic,1)./numel(beh_use_vec);
        %% get the finger print for each type
        behavior_vector(kk,un_el+1) = counts;
        fprintf('amount of time unstructured %s is %f \n',analysisstruct.conditionnames{annotation_choose(kk)},behavior_vector(kk,1));        
    end
    
        behavior_composition_vector{ll} = behavior_vector;
end

hierarchystruct.clustered_behavior = clustered_behavior;
hierarchystruct.annotation_choose = annotation_choose;%([1,2,3,7]);
hierarchystruct.base_annotation = base_annotation;
hierarchystruct.base_annotation_inds = base_annotation_inds;
hierarchystruct.annotation_indexes=annotation_indexes;
hierarchystruct.timescales=timescales;
hierarchystruct.behavior_composition_vector = behavior_composition_vector;
hierarchystruct.normshades_agg = normshades_agg;

end

