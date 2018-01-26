clear all
close all
clc

inputFile = '../data/comb/comb';

pwh_list = readPWHList([inputFile, '.dat']);

figure
hold on
plotPWHList(pwh_list)
axis equal
title('Original')

outputIdentifier = 'test';

figure

filename = [inputFile, '_', outputIdentifier, '_original_close'];
pwh_list = readPWHList([filename, '.dat']);

subplot(2, 2, 1)
hold on
plotPWHList(pwh_list)
axis equal
title('Original Close')

filename = [inputFile, '_', outputIdentifier, '_original_open'];
pwh_list = readPWHList([filename, '.dat']);

subplot(2, 2, 2)
hold on
plotPWHList(pwh_list)
axis equal
title('Original Open')

filename = [inputFile, '_', outputIdentifier, '_original_thin'];
pwh_list = readPWHList([filename, '.dat']);

subplot(2, 2, 3)
hold on
plotPWHList(pwh_list)
axis equal
title('Original Thin')

filename = [inputFile, '_', outputIdentifier, '_original_hole'];
pwh_list = readPWHList([filename, '.dat']);

subplot(2, 2, 4)
hold on
plotPWHList(pwh_list)
axis equal
title('Original Hole')

filename = [inputFile, '_', outputIdentifier];
pwh_list = readPWHList([filename, '.dat']);

figure
hold on
plotPWHList(pwh_list)
axis equal
title('SDFix Output')

figure

filename = [inputFile, '_', outputIdentifier, '_sdfix_close'];
pwh_list = readPWHList([filename, '.dat']);

subplot(2, 2, 1)
hold on
plotPWHList(pwh_list)
axis equal
title('SDFix Close')

filename = [inputFile, '_', outputIdentifier, '_sdfix_open'];
pwh_list = readPWHList([filename, '.dat']);

subplot(2, 2, 2)
hold on
plotPWHList(pwh_list)
axis equal
title('SDFix Open')

filename = [inputFile, '_', outputIdentifier, '_sdfix_thin'];
pwh_list = readPWHList([filename, '.dat']);

subplot(2, 2, 3)
hold on
plotPWHList(pwh_list)
axis equal
title('SDFix Thin')

filename = [inputFile, '_', outputIdentifier, '_sdfix_hole'];
pwh_list = readPWHList([filename, '.dat']);

subplot(2, 2, 4)
hold on
plotPWHList(pwh_list)
axis equal
title('SDFix Hole')

figure

filename = [inputFile, '_', outputIdentifier, '_final_close'];
pwh_list = readPWHList([filename, '.dat']);

subplot(2, 2, 1)
hold on
plotPWHList(pwh_list)
axis equal
title('Final Close')

filename = [inputFile, '_', outputIdentifier, '_final_open'];
pwh_list = readPWHList([filename, '.dat']);

subplot(2, 2, 2)
hold on
plotPWHList(pwh_list)
axis equal
title('Final Open')

filename = [inputFile, '_', outputIdentifier, '_final_thin'];
pwh_list = readPWHList([filename, '.dat']);

subplot(2, 2, 3)
hold on
plotPWHList(pwh_list)
axis equal
title('Final Thin')

filename = [inputFile, '_', outputIdentifier, '_final_hole'];
pwh_list = readPWHList([filename, '.dat']);

subplot(2, 2, 4)
hold on
plotPWHList(pwh_list)
axis equal
title('Final Hole')

