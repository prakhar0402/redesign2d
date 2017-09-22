function plotPWHvertices(pwh, color, size)

if nargin < 2
    color = [0, 0, 0];
    if nargin < 3
        size = 1;
    end
end

plot(pwh.outer_boundary(:, 1), pwh.outer_boundary(:, 2), 'o',  'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerSize', size);
for i = 1 : pwh.num_holes
    plot(pwh.holes{i}(:, 1), pwh.holes{i}(:, 2), 'o',  'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerSize', size);
end