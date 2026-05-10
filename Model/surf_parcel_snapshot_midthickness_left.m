function surf_parcel_snapshot_midthickness_left(ts,stim,lims,conn_ab)
%SURF_SCHAEF_MOVIE_FLAT_LEFT Animate Schaefer data on flat left hemisphere
%   Input: data must be 400 x time

% Load flat map data
surf = gifti('../Model/Data/fs_LR.32k.L.midthickness.surf.gii');
V = double(surf.vertices); F = double(surf.faces); 
clear surf;

% Load parcellation labels per vertex (400 parcels)
V_P_map = dlmread('../Model/Data/fsLR_32k_Schaefer400-lh.txt');

if isempty(lims)
    lims = max(ts, [], 'all')*[0, 1];
end

nTime = size(ts, 2);
nparcels = size(ts, 1);

% If V_P_map(i) == 0, then ts(V_P_map(i), initial_time) should be 0 i
mask = V_P_map > 0;             % Valid indices

Colormap = [linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)'];

% Set up figure and plot once
ax = gca;
hold;
set(gcf, 'Renderer', 'zbuffer', 'Color', 'w'); % Set background color to white
colormap(ax, Colormap);
axis equal;
axis off;
set(ax, 'XColor', 'none', 'YColor', 'none');  % Hide axis ticks
view(2);
result = zeros(size(V_P_map));  % Preallocate output
result(mask) = ts(V_P_map(mask));
p = patch('Vertices', V, ...
          'Faces', F, ...
          'FaceVertexCData', result, ...% One plus cos first col is the idxs
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');
% p.FaceAlpha = 0.5;
view(90, 0);
% uistack(tl, 'top'); % Bring the text to the top layer
clim(lims);
% cb = colorbar; cb.TickLabelInterpreter = 'latex';
% PLot stimulus positions and lrcs on separate axes on top
light('Position',[0 0 0],'Style','infinite')
material([0.95 0 0]);


% Create second axes on top
ax2 = axes('Position', ax.Position, ...
    'Color', 'none', ...
    'XLim', ax.XLim, 'YLim', ax.YLim, 'ZLim', ax.ZLim, ...
    'CameraPosition', ax.CameraPosition, ...
    'CameraTarget', ax.CameraTarget, ...
    'CameraUpVector', ax.CameraUpVector, ...
    'DataAspectRatio', ax.DataAspectRatio, ...
    'Visible', 'off');

% Draw your black cross here, same coordinates as original
axes(ax2); hold on;
for pos = stim
    % Get vertices for parcel 1
    stim_P = find(V_P_map == pos);
    % Compute centroid of parcel 1
    stim_P_centroid = mean(V(stim_P, :), 1);
    % Plot cross
    cross_handle = plot3(stim_P_centroid(1), stim_P_centroid(2), stim_P_centroid(3), ...
        'Color', 'k', 'Marker', '.', 'MarkerSize', 20);
    uistack(cross_handle, 'top');
end

if ~(isempty(conn_ab))
% Find coordinates of parcels 30 and 100
    for i = 1:size(conn_ab, 1)
        % Compute centroid positions
        a_P = find(V_P_map == conn_ab(i, 1));
        a = mean(V(a_P, :), 1);
        b_P = find(V_P_map == conn_ab(i, 2));
        b = mean(V(b_P, :), 1);
        % Compute direction vector
        dir = b - a;
        % Draw an arrow from a to b
        quiver3( ...
        a(1), a(2), a(3), ...
            dir(1), dir(2), dir(3), ...
            0, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5, 'MaxHeadSize', 0);
        drawnow;
    end
end

hold off;