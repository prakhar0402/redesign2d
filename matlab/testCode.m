clear all
close all
clc

pwh_list = readPWHList('../data/exampleOut.dat');

normals = -dlmread('../data/normals.dat');
normal_length = 0.5;

angles = dlmread('../data/angles.dat')*180/pi;

figure
hold on
plotPWHList(pwh_list)
axis equal

count = 0;
for i = 1 : size(pwh_list.pwh{1}.outer_boundary, 1);
    
    v = pwh_list.pwh{1}.outer_boundary(i, :);
    count = count + 1;
    x = [v(1), v(1)+normal_length*normals(count, 1)];
    y = [v(2), v(2)+normal_length*normals(count, 2)];
    plot(x, y, 'LineWidth', 3)
end
for h = 1 : pwh_list.pwh{1}.num_holes
    for i = 1 : size(pwh_list.pwh{1}.holes{h}, 1);

        v = pwh_list.pwh{1}.holes{h}(i, :);
        count = count + 1;
        x = [v(1), v(1)+normal_length*normals(count, 1)];
        y = [v(2), v(2)+normal_length*normals(count, 2)];
        plot(x, y, 'LineWidth', 3)
    end
end

sdf = dlmread('../data/SDF.dat');

figure
hold on
plotPWHwithCdata(pwh_list.pwh{1}, sdf)
axis equal
colorbar