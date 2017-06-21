clear all
close all
clc
inputFile = '../data/example/example';

% pwh_list = readPWHList([inputFile, '.dat']);

outputIdentifier = '0003';
filename = [inputFile, '_', outputIdentifier];

open_pwh_list = readPWHList([filename, '_open.dat']);

figure
hold on
plotPWHList(open_pwh_list)
axis equal

close_pwh_list = readPWHList([filename, '_close.dat']);

figure
hold on
plotPWHList(close_pwh_list)
axis equal
