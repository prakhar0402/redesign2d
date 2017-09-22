function pwh_list = readPWHList(filename)

% PWH List MAT Data Structure
% num_poly: Number of pwhs
% pwh: cell array containing pwhs
% pwh{i}
%       num_holes: number of holes in pwh
%       holes: cell array containing holes
%       outer_boundary: Nx2 array of vertices
%       holes{j}
%               Mx2 array of vertices

fid = fopen(filename);
pwh_list.num_poly = fscanf(fid, '%d', 1);
pwh_list.pwh = cell(pwh_list.num_poly, 1);
for i = 1 : pwh_list.num_poly
    pwh_list.pwh{i}.num_holes = fscanf(fid, '%d', 1);
    pwh_list.pwh{i}.holes = cell(pwh_list.pwh{i}.num_holes, 1);
    num_V = fscanf(fid, '%d', 1);
    pwh_list.pwh{i}.outer_boundary = fscanf(fid, '%f', [2 num_V])';
    for j = 1 : pwh_list.pwh{i}.num_holes
        num_V = fscanf(fid, '%d', 1);
        pwh_list.pwh{i}.holes{j} = fscanf(fid, '%f', [2 num_V])';
    end
end
fclose(fid);