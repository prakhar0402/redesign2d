function plotPWHwithCdata(pwh, cdata, edgecolor, facealpha, holealpha)

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

count = size(pwh.outer_boundary, 1);
patch('Faces', 1:count, 'Vertices', pwh.outer_boundary, 'FaceVertexCData', cdata(1:count), 'FaceColor', 'interp', 'EdgeColor', edgecolor, 'FaceAlpha', facealpha);
for i = 1 : pwh.num_holes
    patch('Faces', 1:size(pwh.holes{i}, 1), 'Vertices', pwh.holes{i}, 'FaceVertexCData', cdata(count+1:count+size(pwh.holes{i}, 1)), 'FaceColor', 'interp', 'EdgeColor', edgecolor, 'FaceAlpha', holealpha);
    count = count + size(pwh.holes{i}, 1);
end