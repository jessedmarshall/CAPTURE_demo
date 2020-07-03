function [fulloutput,indivbouts] = fillannotationgaps(inputarray,numframes,add_edges)

inputarray = unique(sort(inputarray,'ASCEND'));

diff_input = diff(inputarray);
fulloutput = [];
%define borders
diff_input(diff_input>numframes) = 0;
diff_input = [numframes diff_input];
outputstruct = bwconncomp(diff_input);
indivbouts = cell(1,numel(outputstruct.PixelIdxList));

border_amount = ceil(numframes./2);
maxval =max(inputarray);
for mm = 1:numel(outputstruct.PixelIdxList)
    if (add_edges)
        good_inds = (inputarray( outputstruct.PixelIdxList{mm}(1))-border_amount):(inputarray(outputstruct.PixelIdxList{mm}(end))+border_amount);
        good_inds(good_inds<1) = 1;
        good_inds(good_inds>maxval) = maxval;
        
           fulloutput = cat(2,fulloutput,good_inds);
   indivbouts{mm} = good_inds;
   
    else
   fulloutput = cat(2,fulloutput,inputarray( outputstruct.PixelIdxList{mm}(1)):inputarray(outputstruct.PixelIdxList{mm}(end)));
   indivbouts{mm} = inputarray( outputstruct.PixelIdxList{mm}(1)):inputarray(outputstruct.PixelIdxList{mm}(end));
    end
end




end