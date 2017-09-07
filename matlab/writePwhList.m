function writePwhList(fileId, pwh_list)

fprintf(fileId, '%d\n', pwh_list.num_poly);
for i = 1 : pwh_list.num_poly
    writePwh(fileId, pwh_list.pwh{i}, 1);
end