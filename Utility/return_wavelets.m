function [wavelet_coeffs,wtAll,fr]= WaveletCluster(features, frames, opts ,wtpre,frpre)

% DEPENDENCIES:
% addpath(genpath(...\Movement_analysis-master\MotionMapper-master\));
% addpath(genpath(...\Movement_analysis-master\VbGm));
%addpath(genpath('D:\Kevin\repositories2\behavior-mapping\MotionMapper'));
%addpath(genpath('D:\Kevin\Motion_mapping\utilities'));
%addpath('D:\Kevin\Motion_mapping\VbGm');
% clustering script using wavelet space, and GMM / vb method
%
% features is T x d (observations x dimensions)
%  - this assumes d is relatively small (~<50)!!!
%

%% initialize options
if nargin<3 || isempty(opts);
    frames = [];
    
    opts.whiten = 0;
    opts.frameNormalize = 0;
    opts.clustermethod = 'GMM';
    num = size(features,1); % number modes (spectrograms to find) (usually want full dimension)
    % num = pcuse;
    ds = 1; % down sampling
	samprate = 245;
    params = struct;
    params.samplingFreq = samprate/ds;
    params.numPeriods=50; %distinct number of frequencies to use
    params.minF = 0.2; % min freq to analyze
    params.maxF = 30; % max freq to analyze
    
    
    % for simulated data;
    params.samplingFreq = 500;
    
    
    % for gmm model parameters
    pcuse = 20;
    numclusters = 40;
    lambda = 0.1; % regularization
   
    
else
    num=opts.num;
    ds = opts.ds;
    samprate = opts.samprate;
    params = opts.params;
    
    pcuse = opts.pcuse;
    lambda = opts.lambda;
end
 
    subsample = 10; % number of times to run GMM fit on subsample
    numsample = 500; % number of data points to sample

%% crop frames
if ~isempty(frames);
features = features(frames,:);
end
T = size(features,1);
d = size(features,2);

%% whiten the data
if opts.whiten
features = bsxfun(@rdivide,bsxfun(@minus,features,mean(features,2)),std(features,[],2));
end

%% find wavelet transform (dependent on berman code!!!)
%params.omega0 = 5;
%params.samplingFreq = 500;
if nargin == 5
    disp('already using precalculated wavelet');
    wt = wtpre;
    fr = frpre;
else
    disp('finding wavelet transform');
    [wt,fr]=findWavelets(features,num,setRunParameters(params));
end
wt = log(wt);

% cut out some freqs
%wt = wt(:,fr > 4);
%fr = fr(fr > 4);

%figure; 
%ax1=subplot(3,1,[1 2]);imagesc(wt'); %set(ax1,'visible','off');
%title(['omega = ' num2str(params.omega0)]);
%ax2 = subplot(3,1,3); plot(features); set(ax2,'xlim',[0 length(features)]);
%linkaxes([ax1 ax2],'x');

% wt should be T x d*numPeriods

if opts.frameNormalize;

% frame normalize (see Todd paper)
wtAmps = sum(wt,2);
% exclude low energy frames
options = statset('MaxIter',1000);
gmmodel = fitgmdist(wtAmps,5,'Options',options,'replicates',1);
% thresh as min between two largest peaks
range = quantile(wtAmps,[.01 .99]);
range = linspace(range(1)-std(wtAmps),range(2)+std(wtAmps),1000);
gaus_fit = pdf(gmmodel,range');
[pks,loc] = findpeaks(gaus_fit);
[~,threshind] = min(gaus_fit(loc(1):loc(2))); threshind = threshind + loc(1);
% figure;plot(range,gaus_fit);
% hold on;plot(range(loc),pks,'ro');
% plot(range(threshind),gaus_fit(threshind),'go');

thresh = range(threshind);
useind = find(wtAmps > thresh);
keepind = find(wtAmps < thresh);

%figure; hold on;
%histogram(wtAmps(wtAmps > thresh)); histogram(wtAmps(wtAmps < thresh));

% frame normalize
wtNorm = bsxfun(@rdivide,wt,wtAmps);
wtNorm = wtNorm .* max(max(wt))/max(max(wtNorm)); % rescale
% figure;imagesc(wtNorm(:,1:100)');colorbar; title('wt norm');
% figure;imagesc(wt(:,1:100)');colorbar; title('wt raw');
wtAll = wtNorm;
wtAll(keepind,:) = wt(keepind,:);
% figure;imagesc(wtAll(:,1:100)');colorbar; title('wt all');
else
wtAll = wt; % because frame normalization still a bit sketchy
end


%% take PCA of wavelets
disp('finding PCs');
    wtAll(find(isinf(wtAll))) = 0;
     wtAll(find(isnan(wtAll))) = 0;
   
if exist('pca_randomized')
[U,S,V] = pca_randomized(wtAll, min(pcuse,size(wtAll,1)));
coeff = zeros(size(wtAll,2),size(wtAll,2));

score = U;
coeff = S;
% use only if don't have pca_randomized]

explained = diag(S).*diag(S) ./ sum(diag(S).*diag(S));

else
[coeff,score,latent,tsquared,explained,mu] = pca(wtAll);
[U,S,V] = svd(wtAll);
end

feat_pcs = V;
wavelet_coeffs = U;
% 
% figure(667)
% plot3(wavelet_coeffs(1:50:500000,4),wavelet_coeffs(1:50:500000,2),wavelet_coeffs(1:50:500000,3),'+')
% % weight pca modes to prevent noise from overwhelming GMM


% plots wavelet feature space sorted into clusters
% [sorted_labels,ind] = sort(labels,'ascend');
% figure;imagesc(wtAll(labels,:));
% 
% label_ind = find(diff(sorted_labels)>0);
% hold on;
% for j = 1:length(label_ind);
%     line([0 size(wtAll,2)],[label_ind(j) label_ind(j)],'color','r');
% end

end
