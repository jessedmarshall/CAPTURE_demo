function [D,entropies] = findListKLDivergences(data,data2)
%finds the KL-divergences (D) and entropies between all rows in 'data' and
%all rows in 'data2'
    
    % JT: adding eps here to avoid generating -infs, note that there is similar logic in findKLDivergences (infs set to zero in that case)
   % logData = log(data+eps);

   % N = length(data(:,1));
    logData = log(data);
    logData(isinf(logData) | isnan(logData)) = 0;
    
    
    entropies = -sum(data.*logData,2);
    clear logData;  

    logData2 = log(data2+eps);  

    D = - data * logData2';
    
    D = bsxfun(@minus,D,entropies); 
        
    D = D ./ log(2);
