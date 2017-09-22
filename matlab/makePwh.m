clear all
close all
clc

%% case 1
% pwh_list.num_poly = 1;
% pwh_list.pwh = cell(pwh_list.num_poly, 1);
% pwh_list.pwh{1}.num_holes = 2;
% pwh_list.pwh{1}.holes = cell(pwh_list.pwh{1}.num_holes, 1);
% pwh_list.pwh{1}.outer_boundary = csvread('../data/case1/case1_o.csv');
% pwh_list.pwh{1}.holes{1} = csvread('../data/case1/case1_h1.csv');
% pwh_list.pwh{1}.holes{2} = csvread('../data/case1/case1_h2.csv');
% 
% fid = fopen('../data/case1/case1.dat', 'w');
% writePwhList(fid, pwh_list);
% fclose(fid);

%% case 2
% pwh_list.num_poly = 1;
% pwh_list.pwh = cell(pwh_list.num_poly, 1);
% pwh_list.pwh{1}.num_holes = 1;
% pwh_list.pwh{1}.holes = cell(pwh_list.pwh{1}.num_holes, 1);
% pwh_list.pwh{1}.outer_boundary = csvread('../data/case2/case2_o.csv');
% pwh_list.pwh{1}.holes{1} = csvread('../data/case2/case2_h.csv');
% 
% fid = fopen('../data/case2/case2.dat', 'w');
% writePwhList(fid, pwh_list);
% fclose(fid);

%% sample 2
% pwh_list.num_poly = 1;
% pwh_list.pwh = cell(pwh_list.num_poly, 1);
% pwh_list.pwh{1}.num_holes = 1;
% pwh_list.pwh{1}.holes = cell(pwh_list.pwh{1}.num_holes, 1);
% pwh_list.pwh{1}.outer_boundary = csvread('../data/sample2/sample2_o.csv');
% pwh_list.pwh{1}.holes{1} = csvread('../data/sample2/sample2_h.csv');
% 
% fid = fopen('../data/sample2/sample2.dat', 'w');
% writePwhList(fid, pwh_list);
% fclose(fid);

%% case 3
% pwh_list.num_poly = 1;
% pwh_list.pwh = cell(pwh_list.num_poly, 1);
% pwh_list.pwh{1}.num_holes = 4;
% pwh_list.pwh{1}.holes = cell(pwh_list.pwh{1}.num_holes, 1);
% pwh_list.pwh{1}.outer_boundary = csvread('../data/case3/case3_o.csv');
% pwh_list.pwh{1}.holes{1} = csvread('../data/case3/case3_h1.csv');
% pwh_list.pwh{1}.holes{2} = csvread('../data/case3/case3_h2.csv');
% pwh_list.pwh{1}.holes{3} = csvread('../data/case3/case3_h3.csv');
% pwh_list.pwh{1}.holes{4} = csvread('../data/case3/case3_h4.csv');
% 
% fid = fopen('../data/case3/case3.dat', 'w');
% writePwhList(fid, pwh_list);
% fclose(fid);

%% case 4
pwh_list.num_poly = 1;
pwh_list.pwh = cell(pwh_list.num_poly, 1);
pwh_list.pwh{1}.num_holes = 0;
pwh_list.pwh{1}.holes = cell(pwh_list.pwh{1}.num_holes, 1);
pwh_list.pwh{1}.outer_boundary = csvread('../data/case4/case4_o.csv');

fid = fopen('../data/case4/case4.dat', 'w');
writePwhList(fid, pwh_list);
fclose(fid);

%% display
figure
hold on
plotPWHList(pwh_list)
axis equal