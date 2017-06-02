function plotPWH(pwh, facecolor, edgecolor, facealpha, holealpha)

if nargin < 5
    holealpha = 1;
    if nargin < 4
        facealpha = 1;
        if nargin < 3
            edgecolor = [0 0 0];
            if nargin < 2
                facecolor = [1 0 0];
            end
        end
    end
end

patch('Faces', 1:size(pwh.outer_boundary, 1), 'Vertices', pwh.outer_boundary, 'FaceColor',  facecolor, 'EdgeColor', edgecolor, 'FaceAlpha', facealpha);
for i = 1 : pwh.num_holes
    patch('Faces', 1:size(pwh.holes{i}, 1), 'Vertices', pwh.holes{i}, 'FaceColor', [1, 1, 1], 'EdgeColor', edgecolor, 'FaceAlpha', holealpha);
end