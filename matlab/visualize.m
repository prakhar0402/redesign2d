clear all
close all
clc

% inputFile = '../data/example/example';
% inputFile = '../data/cactus_z0/cactus_z0';
% inputFile = '../data/karambit/karambit_2_1';
% inputFile = '../data/sample/sample';
% inputFile = '../data/case4/case4';
% inputFile = '../data/sample2/sample2';
inputFile = '../data/ex4/ex4';

% pwh_list = readPWHList([inputFile, '.dat']);

% outputIdentifier = '0000';
outputIdentifier = 'test';
filename = [inputFile, '_', outputIdentifier];

pwh_list = readPWHList([filename, '.dat']);
fid = fopen([filename, '.out']);

line = fgetl(fid);
fileloc = fgetl(fid);
line = fgetl(fid);
identifier = fgetl(fid);
line = fgetl(fid);
T1 = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
K_S = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
K_D = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
K_SDF = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
K_C = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
VERTEX_MASS = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
TIME_STEP = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
TOTAL_TIME = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
STEPS = fscanf(fid, '%d\n', 1);
line = fgetl(fid);
nMove = fscanf(fid, '%d\n', 1);
line = fgetl(fid);
nV = fscanf(fid, '%d\n', 1);
line = fgetl(fid);
sdf = fscanf(fid, '%f\n', nV);
line = fgetl(fid);
movable_index = fscanf(fid, '%f\n', nV);
line = fgetl(fid);
normals = fscanf(fid, '%f %f\n', [2 nV])';
line = fgetl(fid);
angles = fscanf(fid, '%f %f\n', [2 nV])';
line = fgetl(fid);
interiorAngles = fscanf(fid, '%f\n', nV);
line = fgetl(fid);
max_change = fscanf(fid, '%f\n', STEPS);
line = fgetl(fid);
init_rho = fscanf(fid, '%f\n', 1);
line = fgetl(fid);
sdfix_rho = fscanf(fid, '%f\n', 1);

fclose(fid);

figure
subplot(2, 2, 1)
hold on
plotPWHwithCdata(pwh_list.pwh{1}, sdf)
colorbar
axis equal
xlabel({['SDF: ', inputFile], outputIdentifier}, 'fontweight', 'bold')
set(gca, 'XTick', '', 'YTick', '')
box on

% saveas(gca, [filename, '.jpg'])

% T1 = 6.2;

remain_frac = sum(sdf < T1)/nV;
fprintf('Fraction of remaining vertices below threshold diameter = %f\n', remain_frac);

subplot(2, 2, 2)
plotPWHwithCdata(pwh_list.pwh{1}, 1*(sdf < T1) + 0*(movable_index > 0))
colorbar
axis equal
xlabel({['Thin: ', inputFile], outputIdentifier}, 'fontweight', 'bold')
set(gca, 'XTick', '', 'YTick', '')
box on
hold on



% plotting normals
subplot(2, 2, 3)
hold on
plotPWHList(pwh_list)

normal_length = 0.05;
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
axis equal
xlabel({['Normals: ', inputFile], outputIdentifier}, 'fontweight', 'bold')
set(gca, 'XTick', '', 'YTick', '')
box on
hold on

subplot(2, 2, 4)
plotPWHwithCdata(pwh_list.pwh{1}, interiorAngles*180/pi)
colorbar
axis equal
xlabel({['Angles: ', inputFile], outputIdentifier}, 'fontweight', 'bold')
set(gca, 'XTick', '', 'YTick', '')
box on
hold on

%% Convergence plot
figure
plot(TIME_STEP:TIME_STEP:TIME_STEP*STEPS, max_change)
xlabel('Time (in secs)');
ylabel('Maxima of Change in Vertex Location')

%% Open and Close Plot
open_pwh_list = readPWHList([filename, '_open.dat']);

figure
subplot(1, 2, 1)
hold on
plotPWHList(open_pwh_list)
axis equal
xlabel({['Open: ', inputFile], outputIdentifier}, 'fontweight', 'bold')
set(gca, 'XTick', '', 'YTick', '')
box on
hold on

close_pwh_list = readPWHList([filename, '_close.dat']);

subplot(1, 2, 2)
hold on
plotPWHList(close_pwh_list)
axis equal
xlabel({['Close: ', inputFile], outputIdentifier}, 'fontweight', 'bold')
set(gca, 'XTick', '', 'YTick', '')
box on
hold on