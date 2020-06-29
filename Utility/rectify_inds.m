function array_out = rectify_inds(array_in,max_length)
array_in(array_in<1) = 1;
array_in(array_in>max_length) = max_length;
array_out = array_in;
end