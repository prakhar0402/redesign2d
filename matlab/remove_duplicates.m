close all
clear all
clc

filename = '../data/comb/comb.xlsx';

data = xlsread(filename);
data = [data; data(1, :)];
new_data = data(1, :);

for i = 2 : size(data, 1)
    d = sum((data(i, :) - new_data(end, :)).^2);
    if d > 0
        new_data = [new_data; data(i, :)];
    end
end

new_data = 1000*new_data;

filename = '../data/comb/comb_nd_x1000.xlsx';
xlswrite(filename, new_data)
        