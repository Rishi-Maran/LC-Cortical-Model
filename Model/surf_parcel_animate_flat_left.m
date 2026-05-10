function surf_parcel_animate_flat_left(ts,dt,t0,t1,step,stim,lims,isgif,conn_ab)
%SURF_SCHAEF_MOVIE_FLAT_LEFT Animate Schaefer data on flat left hemisphere
%   Input: data must be 400 x time

% Load flat map data
surf = gifti('../Model/Data/fs_LR.32k.L.flat.surf.gii');
V = double(surf.vertices); F = double(surf.faces); 
clear surf;

% Load parcellation labels per vertex (400 parcels)
V_P_map = dlmread('../Model/Data/fsLR_32k_Schaefer400-lh.txt');

initial_time = max(1, floor(t0 / dt));
final_time = min(size(ts, 2), ceil(t1 / dt));

ts = ts(:, initial_time:final_time);

if nargin < 7 || isempty(lims)
    lims = max(ts, [], 'all')*[-1, 1];
end

nTime = size(ts, 2);
nparcels = size(ts, 1);

% If V_P_map(i) == 0, then ts(V_P_map(i), initial_time) should be 0 i
mask = V_P_map > 0;             % Valid indices

% ts_v = zeros()

% Build lookup matrix
% idx = zeros(1+nparcels, 1 + nTime);
% idx(2:end, 1) = 1:nparcels;
% idx(2:end, 2:end) = ts;

% Precompute index mapping
% original = map;
% [~, index_L] = ismember(original, idx(:, 1));

Colormap = [linspace(0, 1, 256)' linspace(0, 1, 256)' linspace(1, 1, 256)';
    linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)'];

% Set up figure and plot once
f = figure; f.Position = [100, 100, 800, 600];
hold;
set(gcf, 'Renderer', 'zbuffer', 'Color', 'w'); % Set background color to white
colormap(Colormap);
axis equal;
set(gca, 'XColor', 'none', 'YColor', 'none');  % Hide axis ticks
view(2);
result = zeros(size(V_P_map));  % Preallocate output
result(mask) = ts(V_P_map(mask), initial_time);
p = patch('Vertices', V, ...
          'Faces', F, ...
          'FaceVertexCData', result, ...% One plus cos first col is the idxs
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');
% title(gca, Title,'FontSize',20);
tl = text(0.5, 0, sprintf('Time = %.2f s', t0), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 15, 'FontWeight', 'normal', 'units', 'normalized');
% uistack(tl, 'top'); % Bring the text to the top layer
clim(lims);
colorbar;

% --- Detect boundary edges ---
E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];           % all edges
E_sorted = sort(E, 2);                              % ensure consistent ordering
[uniqueEdges, ~, ic] = unique(E_sorted, 'rows');    % unique edge list
counts = histcounts(ic, 1:(max(ic)+1));             % count how many times each edge appears
boundaryEdges = uniqueEdges(counts == 1, :);        % keep edges that appear only once

% --- Plot boundary edges ---
hold on
for i = 1:size(boundaryEdges, 1)
    v1 = V(boundaryEdges(i,1), :);
    v2 = V(boundaryEdges(i,2), :);
    plot3([v1(1), v2(1)], [v1(2), v2(2)], [v1(3), v2(3)], ...
          'k-', 'LineWidth', 1);  % black line
end

for pos = stim
    % Get vertices for parcel 1
    stim_P = find(V_P_map == pos);
    % Compute centroid of parcel 1
    stim_P_centroid = mean(V(stim_P, :), 1);
    % Plot an X marker at the centroid
    xMarker = plot3(stim_P_centroid(1), ...
        stim_P_centroid(2), stim_P_centroid(3), ...
        'kx', 'MarkerSize', 10, 'LineWidth', 2);
end

if ~(isempty(conn_ab))
% Find coordinates of parcels 30 and 100
    for i = 1:size(conn_ab, 1)
        % Compute centroid positions
        a = mean(V(V_P_map == conn_ab(i, 1), :), 1);
        b = mean(V(V_P_map == conn_ab(i, 2), :), 1);
        % Compute direction vector
        dir = b - a;
        % Draw an arrow from a to b
        quiver3( ...
        a(1), a(2), a(3), ...
            dir(1), dir(2), dir(3), ...
            0, 'Color', 'k', 'LineWidth', 2, 'MaxHeadSize', 0.1);
        drawnow;
    end
end



% pause;

if isgif

    filename = 'myAnimation.gif';

    for t = initial_time:step:final_time
        result(mask) = ts(V_P_map(mask), t);
        p.FaceVertexCData = result;
        set(tl, 'String', sprintf('Time = %.2f s', t * dt), 'FontSize', 15);
        drawnow;
        frame = getframe(gca);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if t == initial_time
            % First frame: create the GIF and set loop count
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.);
        else
            % Subsequent frames: append to the GIF
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.);
        end    
    end
else
    % Loop through time
    % while(1)
        for t = initial_time:step:final_time
            result(mask) = ts(V_P_map(mask), t);
            p.FaceVertexCData = result;
            set(tl, 'String', sprintf('Time = %.2f s', t * dt), 'FontSize', 15);
            drawnow;
        end
    % end
end