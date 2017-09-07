function writePwh(fileId, pwh, isInList)

if (~isInList)
    fprintf(fileId, '1\n');
end
fprintf(fileId, '%d\n', pwh.num_holes);
writePg(fileId, pwh.outer_boundary, 1);
for i = 1 : pwh.num_holes
    writePg(fileId, pwh.holes{i}, 1);
end