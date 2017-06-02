function plotPWHList(pwh_list, facecolor, edgecolor, facealpha, holealpha)

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

tf = ishold;
hold on
for i = 1 : pwh_list.num_poly
    plotPWH(pwh_list.pwh{i}, facecolor, edgecolor, facealpha, holealpha);
end

if ~tf
    hold off
end
