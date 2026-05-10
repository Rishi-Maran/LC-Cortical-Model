clear; clc;
% Load the left hemisphere mid-thickness surface (GIfTI)
surf = gifti('../Model/Data/fs_LR.32k.L.midthickness.surf.gii');

% Load parcellation labels per vertex (from NSBLab BrainEigenmodes github repository)
parc = dlmread('../Model/Data/fsLR_32k_Schaefer400-lh.txt');

% Load medial wall
% (https://github.com/Washington-University/HCPpipelines/blob/master/global/templates/standard_mesh_atlases/L.atlasroi.32k_fs_LR.shape.gii)
% (https://www.mail-archive.com/hcp-users@humanconnectome.org/msg02150.html)
medial = gifti('../Model/Data/L.atlasroi.32k_fs_LR.shape.gii');
%
is_medial = (medial.cdata == 1);
v_data = surf.vertices(is_medial);
v_parc_map = parc(is_medial);

% check elements of v_parc_map are between 0 and integer
num_p = max(v_parc_map);
if ~all(mod(v_parc_map, 1) == 0) || ~all(v_parc_map >= 0)
    error('v_parc_map must be integers betwee 0 and sth')
end


% Get dimensions of HCP connectome data
num_v = length(v_data);

% Load high resolution connectome
mf = matfile('../Model/Data/S255-v73.mat'); % this is the original file S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat converted to version -v7.3
sz = size(mf, 'avgSC_L');  % 'A' is the variable name
if ~all(sz == num_v)
    error('wrong dimensions');
end

% Compute coarse grained connectome under Schaefer-400 parcellation
sc_p = zeros(num_p, num_v);

chunkSize = 500;

for rowStart = 1:chunkSize:num_v
    rowEnd = min(rowStart + chunkSize - 1, num_v);
    rows = rowStart:rowEnd;

    % Load chunk of A
    sc_chunk = mf.avgSC_L(rows, :);  % chunkSize x N
    
    for idx = 1:num_p
        sc_p(idx, :) = sc_p(idx, :) + sum(sc_chunk(v_parc_map(rows) == idx, :), 1);
    end
    
end
%
sc_p0 = sc_p;
sc_p = zeros(num_p, num_p);


for idx = 1:num_p

    sc_p(:, idx) = sum(sc_p0(:, v_parc_map(1:num_v) == idx), 2);

end

%% save parcellated struct_mat

% uncomment below if you want to save parcellated connectome
save('../Model/Data/struct_mat_p.mat', 'sc_p');


%% Load created parcellated connectome

clear; clc;
% Load the left hemisphere mid-thickness surface (GIfTI)
v_data = gifti('../Model/Data/fs_LR.32k.L.midthickness.surf.gii').vertices;

% Load parcellation labels per vertex (400 parcels)
v_p_map = dlmread('../Model/Data/fsLR_32k_Schaefer400-lh.txt');

% Load structural connectivity matrix (optional, just to confirm size)
load('../Model/Data/struct_mat_p.mat');


%% Check data dimensions

size(v_data) % should be ~ [32492 3]
length(v_p_map) % should be same number of vertices

num_v = size(v_data);
num_p = size(sc_p, 1);  % 400

size(sc_p) % should be [200 200]

%% Compute centroids of each parcel (only non-zero parcel labels)

p_labels = unique(v_p_map);
p_labels(p_labels == 0) = []; % remove zero if present

centroids = zeros(num_p, 3);

for i = 1:num_p
    % Find vertex indices of parcel i
    verts_idx = find(v_p_map == p_labels(i));
    
    % Extract vertex coordinates
    verts_xyz = v_data(verts_idx, :);
    
    % Compute mean xyz coordinate (centroid)
    centroids(i, :) = mean(verts_xyz, 1);
end

% Compute pairwise Euclidean distance matrix between centroids

% Use pdist2 for efficient pairwise Euclidean distances
dist_p = pdist2(centroids, centroids); % [400 x 400]

%% Save the distance matrix to a file

% uncomment this if you want to save the distance matrix
save('../Model/Data/schaefer400_left_euclidean_dist.mat', 'dist_p');

%% Fit exponential function

% Create logical mask for strictly upper triangle (excluding diagonal)
mask_upper = triu(true(size(dist_p)), 1);

% Extract upper triangle elements into vectors
fit_x = dist_p(mask_upper);
fit_y = sc_p(mask_upper);

% Exp Model with parameters A and lambda
exp_model = @(params, xdata) params(1) * exp(-params(2) * xdata);

% Initial guess for [A, lambda]
init_param = [0 0];

% Perform least squares fit
options = optimset('Display','off');
params_fit = lsqcurvefit(exp_model, init_param, fit_x, fit_y, [0 0], [], options);

% Extract fitted A and lambda
A_fit = params_fit(1);
lambda_fit = params_fit(2);

%% Save exponential fit parameters

% uncomment below if you want to save it
A = A_fit; lambda = lambda_fit;
save('../Model/Data/fitting.mat', 'A', 'lambda');
clear A lambda;

%% Identify outliers as FNPs

% set z_threshold
z_max = 3;
d_min = 0;

% Residuals in log space
log_residuals = log(fit_y) - (log(A_fit) - lambda_fit*fit_x);

% Std deviation of residuals in log space
log_res_std = std(log_residuals(log_residuals > -Inf));

%%
% Identify far points
is_far = fit_x > d_min;

% Identify strong SC points relative to multiplicative threshold
is_strong_sc = log_residuals > z_max * log_res_std;

% Selected outliers
idx_outliers = find(is_far & is_strong_sc);

linearIndices = find(mask_upper);
linearIndices_outliers = linearIndices(idx_outliers);

[i_pairs, j_pairs] = ind2sub(size(dist_p), linearIndices_outliers);
log_predicted = zeros(length(i_pairs), 1); 
log_empirical = zeros(length(i_pairs), 1); 
log_empirical_norm = zeros(length(i_pairs), 2);
for k = 1:length(i_pairs)
    log_predicted(k) = exp_model(params_fit, dist_p(i_pairs(k), j_pairs(k)));
    log_empirical(k) = sc_p(i_pairs(k), j_pairs(k));
    log_empirical_norm(k, 1) = log_empirical(k)/sum(sc_p(i_pairs(k), :));
    log_empirical_norm(k, 2) = log_empirical(k)/sum(sc_p(:, j_pairs(k)));
end

% Optional: display number of such pairs
fprintf('Number of outlier pairs: %d\n', length(i_pairs));

for k = 1:length(i_pairs)
    fprintf('Pair (%d, %d): Distance = %.2f mm, SC = %.4e, Normalized SC = %.4e, %.4e, Predicted = %.4e\n', ...
        i_pairs(k), j_pairs(k), dist_p(i_pairs(k), j_pairs(k)), ...
        log_empirical(k), log_empirical_norm(k, 1), log_empirical_norm(k, 2), log_predicted(k));
end

%% Plot with outliers

%
% Plot result
fig = figure; fig.Position = [100, 50, 300, 400];
scatter(fit_x, fit_y, 1, 'k', 'filled'); hold on;
plot(sort(fit_x), exp_model(params_fit, sort(fit_x)), 'b', 'LineWidth', 2);
text(0.8, 0.05, sprintf('$\\mathrm{SC} \\propto \\exp(- %.2f \\, d)$', lambda_fit), 'Color', 'b', 'Interpreter', 'latex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Units', 'normalized');
scatter(fit_x(idx_outliers), fit_y(idx_outliers), 10, [0, 0.5, 0], 'filled')
xlabel('Euclidean Distance ($d$) (mm)', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('Structural Connectivity (SC)', 'Interpreter', 'latex', 'FontSize', 10);
legend('Data', 'Fitted curve', 'FNPs', 'Interpreter', 'latex', 'Location', 'southoutside');
set(gca, 'YScale', 'log');
set(gca, 'TickLabelInterpreter', 'latex');

%%
% save image
exportgraphics(gcf, 'S1_sc_ed.tiff', 'Resolution', 600);
close(fig);

%% Save outliers

outliers = [i_pairs, j_pairs, log_empirical, log_predicted];
save('../Model/Data/outliers.mat', 'outliers');


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
% result(mask) = ts(V_P_map(mask));
p = patch('Vertices', V, ...
          'Faces', F, ...
          'FaceVertexCData', result, ...% One plus cos first col is the idxs
          'FaceColor', [0.95, 0.95, 0.95], ...
          'EdgeColor', 'none', ...
          'FaceAlpha', 0.01, ...
          'EdgeAlpha', 0.01);
p.FaceAlpha = 0.5;
view(-90, 0);
% Plot stimulus positions and lrcs on separate axes on top
light('Style','infinite')
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

load('../Model/Data/outliers.mat');
N_lrcs = 20;
rng(1, 'twister');
perm = randperm(size(outliers, 1), N_lrcs);
conn_ab = ...
    [outliers(perm, [1 2])];

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
        q = quiver3( ...
        a(1), a(2), a(3), ...
            dir(1), dir(2), dir(3), ...
            0, 'Color', [0, 0.5, 0], 'LineWidth', 1, 'MaxHeadSize', 0);
        uistack(q, 'bottom');
    end
end

hold off;

%%

% Optional: Save as SVG
% Uncomment the following lines if you want to save as SVG
% exportgraphics(gcf, 'S1_mesh.tiff', 'Resolution', 1000);
