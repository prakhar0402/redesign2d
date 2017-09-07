% figure for paper
clear all
close all
clc

%%%%
inputFile = '../data/sample/sample';
outputIdentifier = '0000';
filename = [inputFile, '_', outputIdentifier];

% raw outline
pwh_list = readPWHList([inputFile, '.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

% plot vertices
plotPWHListVertices(pwh_list, [0, 0, 1], 2)

figure
% resampled outline
pwh_list = readPWHList([filename, '.dat']);
plotPWHList(pwh_list, [0.8, 0.8, 0.8])
axis equal
axis off
set(gcf,'Color',[1, 1, 1]);

% plot vertices
plotPWHListVertices(pwh_list, [0, 0, 1], 2)