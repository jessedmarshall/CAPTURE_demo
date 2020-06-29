function [yData,betas,P,errors] = run_tSne(data,parameters,useEuclideanDistance)
%run_tSne runs the t-SNE algorithm on an array of normalized wavelet amplitudes
%
%   Input variables:
%
%       data -> Nxd array of wavelet amplitudes (will normalize if
%                   unnormalized) containing N data points
%       parameters -> struct containing non-default choices for parameters
%
%
%   Output variables:
%
%       yData -> N x parameters.num_tsne_dim array of embedding results
%       betas -> Nx1 array of local region size parameters
%       P -> full space transition matrix
%       errors -> P.*log2(P./Q) as a function of t-SNE iteration
%
%
% (C) Gordon J. Berman, 2014
%     Princeton University
    
    % JT: adding euclidean distance
    if ~exist('useEuclideanDistance','var')
        useEuclideanDistance=false;
    end

    
    if nargin < 2
        parameters = [];
    end
    
    parameters = setRunParameters(parameters);

    
    vals = sum(data,2);
    if max(vals) > 1 || min(vals) < 1
        data = bsxfun(@rdivide,data,vals);
    end
    
    fprintf(1,'Finding Distances\n');
    % JT: adding euclidean distance
    if useEuclideanDistance
        D = squareform(pdist(data));
    else
      D =   squareform(pdist(data,'squaredeuclidean'));
          %      D = find_weighted_distance(data);

       % D = findKLDivergences(data);
    end
    
    
    fprintf(1,'Computing t-SNE\n');
    ydata = fast_tsne(D
    [yData,betas,P,errors] = tsne_d(D,parameters);
