% figures for paper
clear all
close all
clc


%%%%
inputFile = '../data/ex1/ex1';
% inputFile = '../data/ex2/ex2';
% inputFile = '../data/ex3/ex3';
% inputFile = '../data/ex4/ex4';
outputIdentifier = 'test';
filename = [inputFile, '_', outputIdentifier];

figure
% raw
pwh_list = readPWHList([inputFile, '.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

figure
% raw
pwh_list = readPWHList([inputFile, '.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);
% plot vertices
plotPWHListVertices(pwh_list, [0, 0, 1], 5)
% xlim([0.2, 0.4])
% ylim([-0.08, 0.2])

figure
% raw
pwh_list = readPWHList([filename, '.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);
% plot vertices
plotPWHListVertices(pwh_list, [0, 0, 1], 5)
% xlim([0.2, 0.4])
% ylim([-0.08, 0.2])


% read out file
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
% plotting normals
hold on
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
plotPWHListVertices(pwh_list, [0, 0, 1], 5)

normal_length = 0.5;
count = 0;
x = pwh_list.pwh{1}.outer_boundary(:, 1);
y = pwh_list.pwh{1}.outer_boundary(:, 2);
u = normals(:, 1);
v = normals(:, 2);
for h = 1 : pwh_list.pwh{1}.num_holes
    x = [x; pwh_list.pwh{1}.holes{h}(:, 1)];
    y = [y; pwh_list.pwh{1}.holes{h}(:, 2)];
end
quiver(x, y, u, v, normal_length)
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);
% xlim([0.2, 0.4])
% ylim([-0.08, 0.2])



figure
% raw
pwh_list = readPWHList([filename, '.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);
% xlim([0.2, 0.4])
% ylim([-0.08, 0.2])


figure
%SDF
plotPWHwithCdata(pwh_list.pwh{1}, sdf)
colorbar
% caxis([0, 1.2])
axis equal
axis off
set(gca,'FontSize', 12)
set(gcf,'Color',[1, 1, 1]);

% thin color
figure
plotPWHwithCdata(pwh_list.pwh{1}, 1*(sdf < T1))
colorbar
axis equal
axis off
set(gca,'FontSize', 12)
set(gcf,'Color',[1, 1, 1]);

% thin
figure
pwh_list = readPWHList([inputFile, '.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);
figure
pwh_list = readPWHList([inputFile, '.dat']);
plotPWHList(pwh_list, 'r')
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

%% Original Slice

% open
figure
pwh_list = readPWHList([filename, '_original_open.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

% close
figure
pwh_list = readPWHList([filename, '_original_close.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

%raw
figure
pwh_list = readPWHList([inputFile, '.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

% slice - open
% figure
hold on
pwh_list = readPWHList([filename, '_original_thin.dat']);
plotPWHList(pwh_list, 'r', 'none')
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

% close - slice
% figure
hold on
pwh_list = readPWHList([filename, '_original_hole.dat']);
plotPWHList(pwh_list, 'm', 'none')
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

% Convergence plot
figure
plot(TIME_STEP:TIME_STEP:TIME_STEP*STEPS, max_change, 'LineWidth', 2)
xlabel('Time Steps', 'FontWeight', 'bold', 'FontSize', 20);
ylabel({'Maxima of','Vertex Displacements'}, 'FontWeight', 'bold', 'FontSize', 20)
axis square



%% Final slice
% open
figure
pwh_list = readPWHList([filename, '_sdfix_open.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

% close
figure
pwh_list = readPWHList([filename, '_sdfix_close.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);


%sdfix
figure
pwh_list = readPWHList([filename, '.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

% final - open
% figure
hold on
pwh_list = readPWHList([filename, '_sdfix_thin.dat']);
plotPWHList(pwh_list, 'r', 'none')
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

% close - final
% figure
hold on
pwh_list = readPWHList([filename, '_sdfix_hole.dat']);
plotPWHList(pwh_list, 'm', 'none')
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

% Convergence plot
figure
plot(TIME_STEP:TIME_STEP:TIME_STEP*STEPS, max_change, 'LineWidth', 2)
xlabel('Time Steps', 'FontWeight', 'bold', 'FontSize', 20);
ylabel({'Maxima of','Vertex Displacements'}, 'FontWeight', 'bold', 'FontSize', 20)
axis square

%% Final slice
% open
figure
pwh_list = readPWHList([filename, '_final_open.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

% close
figure
pwh_list = readPWHList([filename, '_final_close.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);


%final
figure
pwh_list = readPWHList([filename, '_final_close.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

% final - open
% figure
hold on
pwh_list = readPWHList([filename, '_final_thin.dat']);
plotPWHList(pwh_list, 'r', 'none')
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

% close - final
% figure
hold on
pwh_list = readPWHList([filename, '_final_hole.dat']);
plotPWHList(pwh_list, 'm', 'none')
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);