%% Draw mesh of cortex in blue, with green LRCs

% Load flat map data
surf = gifti('../Model/Data/fs_LR.32k.L.midthickness.surf.gii');
V = double(surf.vertices); F = double(surf.faces); 
clear surf;

% Load parcellation labels per vertex (400 parcels)
V_P_map = dlmread('../Model/Data/fsLR_32k_Schaefer400-lh.txt');

% If V_P_map(i) == 0, then ts(V_P_map(i), initial_time) should be 0 i
mask = V_P_map > 0;             % Valid indices

% Set up figure and plot once
f = figure; f.Position = [100, 100, 800, 600];
ax = gca;
hold;
set(gcf, 'Renderer', 'zbuffer', 'Color', 'w'); % Set background color to white
% colormap(ax, Colormap);
axis equal;
axis off;
set(ax, 'XColor', 'none', 'YColor', 'none');  % Hide axis ticks
view(2);
result = zeros(size(V_P_map));  % Preallocate output
p = patch('Vertices', V, ...
          'Faces', F, ...
          'FaceVertexCData', result, ...% One plus cos first col is the idxs
          'FaceColor', [0.5, 0.5, 1], ...
          'EdgeColor', [0.5, 0.5, 1] - 0.5, ...
          'FaceAlpha', 0.01, ...
          'EdgeAlpha', 0.01);
p.FaceAlpha = 0.5;
view(-90, 0);

% Draw thalamus as a sphere inside the cortical mesh
r = 10; % radius of sphere
% Create sphere mesh (21x21 grid by default)
[X, Y, Z] = sphere(50);    % finer resolution with 50 instead of 20
X = r * X;  
Y = r * Y;  
Z = r * Z;  
s = surf2patch(X,Y,Z);
patch('Faces',s.faces,'Vertices',s.vertices, ...
      'FaceColor',[0.2 0.6 1],'EdgeColor',[0.2 0.6 1]*0.7,'FaceAlpha',1);

% Optional: Save as SVG
% Uncomment the following lines if you want to save as SVG
% exportgraphics(gcf, 'mesh_schematic.tiff', 'Resolution', 1000);


%% GEOMETRIC_NETWORK_PATTERN Creates a circular network with distance-based connections
% The strength of connections decreases with distance between nodes

% Parameters - Easily adjustable
n_nodes = 16;           % Number of nodes around the circle
radius = 1;             % Radius of the circle
center = [0, 0];        % Center of the circle
include_center = true;  % Whether to include a central node
max_connection_distance = 1; % Maximum distance for connections (as fraction of diameter)
strength_decay_rate = 2;


% Visual parameters
node_size = 50;         % Size of nodes
node_color = [0.9, 0.9, 0.9];
centre_node_color = [0.3, 0.3, 0.3];
line_color = [0.6, 0.6, 1]; % Red color [R G B]
centre_line_color = [0.5, 0.5, 0.5]; % Red color [R G B]
highlight_color = [0, 0, 1];
background_color = 'white';

% Connection strength parameters
min_line_width = 1.5;   % Minimum line width
max_line_width = 2;     % Maximum line width
min_alpha = 0;        % Minimum transparency
max_alpha = 1;        % Maximum transparency

% Set non-local projections

rng(1, 'twister');
n_lrcs = 3;
min_lrc_length = 0.2; % proportion of number of nodes
lrc_lengths = round(min_lrc_length*n_nodes) + randi(round(n_nodes/2 - min_lrc_length*n_nodes), n_lrcs, 1);
conn_ab = zeros(n_lrcs, 2);
conn_ab(:, 1) = randi(n_nodes, n_lrcs, 1);
conn_ab(:, 2) = mod(conn_ab(:, 1) + lrc_lengths - 1, n_nodes) + 1;

% Generate node positions
% Nodes around the circle
theta = linspace(0, 2*pi, n_nodes + 1);
theta(end) = []; % Remove the duplicate point

node_positions = zeros(n_nodes, 2);
for i = 1:n_nodes
    node_positions(i, :) = center + radius * [cos(theta(i)), sin(theta(i))];
end

% Add center node if requested
if include_center
    node_positions = [node_positions; center];
    total_nodes = n_nodes + 1;
else
    total_nodes = n_nodes;
end

% Calculate distances and create connections
distances = zeros(total_nodes, total_nodes);
for i = 1:total_nodes
    for j = 1:total_nodes
        distances(i, j) = norm(node_positions(i, :) - node_positions(j, :));
    end
end

% Normalize distances (0 to 1)
max_distance = max(distances(:));
normalized_distances = distances / max_distance;

% Create the plot
figure('Color', background_color);
hold on;
axis equal;
axis off;

% Draw long distance connections
for i = 1:n_lrcs
    l = line([node_positions(conn_ab(i,1), 1), node_positions(conn_ab(i,2), 1)], ...
         [node_positions(conn_ab(i,1), 2), node_positions(conn_ab(i,2), 2)], ...
         'Color', [line_color, 0.5], ...
         'LineWidth', max_line_width);
    if i == 1
        l.Color = [highlight_color, 1];
    end
end

% Draw connections
for i = 1:total_nodes
    for j = i+1:total_nodes  % Only draw each connection once
        dist = normalized_distances(i, j);
        
        % Only draw connection if within maximum distance
        if dist <= max_connection_distance && dist > 0
            % Calculate line properties based on distance
            % Closer nodes get thicker, more opaque lines
            strength = exp( - (strength_decay_rate*(dist / max_connection_distance))^2);
            
            % line_width = min_line_width + (max_line_width - min_line_width) * strength;
            line_width = max_line_width;
            alpha = min_alpha + (max_alpha - min_alpha) * strength;
            
            % Draw the line
            l = line([node_positions(i, 1), node_positions(j, 1)], ...
                 [node_positions(i, 2), node_positions(j, 2)], ...
                 'Color', [line_color, alpha], ...
                 'LineWidth', line_width);
            if i == total_nodes || j == total_nodes
                l.Color = [centre_line_color, 0.2];
                l.LineStyle = ':';
                l.LineWidth = max_line_width*1.25;
                uistack(l, 'bottom');
            elseif (i == 1 && j == 2)
                l.Color = [highlight_color];
            end
        end
    end
end

% Draw nodes
scatter(node_positions(:, 1), node_positions(:, 2), node_size, ...
        node_color, 'filled', 'MarkerEdgeColor', 'black');

if include_center
    scatter(center(1), center(2), node_size, centre_node_color, 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 2);
end

% Finalize plot
% title('Geometric Network Pattern', 'FontSize', 14, 'FontWeight', 'bold');
xlim([-1.2, 1.2]);
ylim([-1.2, 1.2]);

%
% Optional: Save as SVG
% Uncomment the following lines if you want to save as SVG
% saveas(gcf, 'ring_schematic.svg', 'svg');


%% Draw sigmoid curve with and without modulation

clear; clc; 
loadparam;

fig = figure; hold;
fig.Position = [100, 100, 450, 150];

v0 = 200;
  
ax = gca;
plot(1:1000, sigmoid(1:1000, 400, 200, 1), 'Color', 'k');
xlabel('Mean Membrane Potential ($V$)', 'Interpreter', 'latex');
ylabel('Mean Firing Rate ($Q$)', 'Interpreter', 'latex');
xticks([]);
yticks([0, 1]);
yticklabels({'0', '$Q_{\mathrm{max}}$'});
ax.TickLabelInterpreter = 'latex';
plot(1:1000, sigmoid(1:1000, 300, (v0 - 300)/(v0 - 400) * 200, 1), 'Color', 'k', 'LineStyle', '--');
scatter(v0, sigmoid(v0, 400, 200, 1), 'k', 'filled');
scatter(400, sigmoid(400, 400, 200, 1), 'k');
scatter(400, sigmoid(400, 300, (v0 - 300)/(v0 - 400) * 200, 1), 'k');

% Optional: Save as SVG
% Uncomment the following lines if you want to save as SVG
% saveas(gcf, '2_modulation.svg', 'svg');
% close(fig);