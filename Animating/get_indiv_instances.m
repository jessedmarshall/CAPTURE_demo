function [valsout] = get_indiv_instances(frameindsin)
vals = zeros(1,max(frameindsin));
vals(frameindsin) = 1;
aa = bwconncomp(vals);
valsout = aa.PixelIdxList;
end