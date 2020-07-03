function annot_out = reorder_annotation_vec(annot_in,reordering)


annot_out = zeros(1,numel(annot_in));

for k =1:max(max(reordering),max(annot_in))
    annot_out(find(annot_in == reordering(k))) = k;
end

end