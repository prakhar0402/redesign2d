function writePg(fileId, vertices, isInPwh)

% vertices should be CCW (CW for holes) ordered list of vertex coordinates

if (~isInPwh)
    fprintf(fileId, '1\n0\n');
end
fprintf(fileId, '%d\n', size(vertices, 1));
for i = 1 : size(vertices, 1)
    fprintf(fileId, '%f %f\n', vertices(i, 1), vertices(i, 2));
end