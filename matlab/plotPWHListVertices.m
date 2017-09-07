function plotPWHListVertices(pwh_list, color, size)

if nargin < 2
    color = [0, 0, 0];
    if nargin < 3
        size = 1;
    end
end

tf = ishold;
hold on
for i = 1 : pwh_list.num_poly
    plotPWHvertices(pwh_list.pwh{i}, color, size);
end

if ~tf
    hold off
end
