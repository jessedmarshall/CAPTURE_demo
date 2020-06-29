function [trainingSetData,trainingSetAmps] = ...
            runEmbeddingSubSampling(flyNames,parameters)
%runEmbeddingSubSampling generates a training set given a set of .mat files
%
%   Input variables:
%
%       projectionDirectory -> directory path containing .mat projection 
%                               files.  Each of these files should contain
%                               an N x pcaModes variable, 'projections'
%       parameters -> struct containing non-default choices for parameters
%
%
%   Output variables:
%
%       trainingSetData -> normalized wavelet training set 
%                           (N x (pcaModes*numPeriods) )
%       trainingSetAmps -> Nx1 array of training set wavelet amplitudes
%       projectionFiles -> cell array of files in 'projectionDirectory'
%
%
% (C) Gordon J. Berman, 2014
%     Princeton University
    
    
    if nargin < 2
        parameters = [];
    end
    parameters = setRunParameters(parameters);
    
        
    
    N = parameters.trainingSetSize;
    L = length(flyNames);
    numPerDataSet = round(N/L);
    numModes = parameters.pcaModes;
    numPeriods = parameters.numPeriods;
    
    % First gather training set contributions for each file in parallel
    trainingSetDataPerFile={};
    trainingSetAmpsPerFile={};
    parfor i=1:L
        % JT: loading fly data by name instead of by path 
        [yData,signalData,signalAmps,~] = ...
                file_embeddingSubSampling(flyNames{i},parameters);
            
        [trainingSetDataPerFile{i},trainingSetAmpsPerFile{i}] = ...
            findTemplatesFromData(signalData,yData,signalAmps,...
                                numPerDataSet,parameters);
            
            
    end

    trainingSetData = zeros(numPerDataSet*L,numModes*numPeriods);
    trainingSetAmps = zeros(numPerDataSet*L,1);
    
    % Now combine contributions from each file
    for i=1:L
        currentIdx = (1:numPerDataSet) + (i-1)*numPerDataSet;
        trainingSetData(currentIdx,:)=trainingSetDataPerFile{i};
        trainingSetAmps(currentIdx)=trainingSetAmpsPerFile{i};
    end
    
