clear all
close all
clc

filename = '../data/open';
pwh_list = readPWHList([filename, '.dat']);

figure
hold on
plotPWHList(pwh_list)
axis equal
