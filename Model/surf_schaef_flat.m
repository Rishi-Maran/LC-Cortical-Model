function surf_schaef_flat(data, conn_ab)
%SURF_SCHAEF       Plots results onto Schaefer Atlas
%
% Note: data must be 400x1

% Load flat surface data
v_data = gifti('../Model/Data/fs_LR.32k.L.midthickness.surf.gii').vertices;

% Load parcellation labels per vertex (400 parcels)
v_p_map = dlmread('../Model/Data/fsLR_32k_Schaefer400-lh.txt');

% identity matrix
idx = zeros(401,2);
idx(2:401,1) = 1:1:400;
idx(2:401,2) = data;

limits = [min(data), max(data)];

% Check for any NaN or Inf values
if limits(1) == limits(2)
    limits = [-1, 1];  % Set a reasonable default
end


% left
original = left.cdata;
[~,index_net] = ismember(original,idx(:,1)); % original(i) occurs in the index_net(i)th position in idx
map_net = idx(:,2);
left.cdata = map_net(index_net);

% right
original = right.cdata;
[~,index_net] = ismember(original,idx(:,1));
map_net = idx(:,2);
right.cdata = map_net(index_net);

% Create figure with subplots
f = figure;
f.Position = [100, 100, 1100, 500];
hold;

% % Left hemisphere subplot
% figure;
% % subplot(1, 2, 1);
% plot(struc_L, left);  % This plots the left hemisphere
% caxis(limits);  % Apply color axis limits
% colorbar;       % Add color bar
% 
% % Right hemisphere subplot
% figure;
% % subplot(1, 2, 2);
% plot(struc_R, right);  % This plots the right hemisphere
% caxis(limits);  % Apply color axis limits
% colorbar;       % Add color bar

% Extract surface data
vertices_L = double(struc_L.vertices);
faces_L = double(struc_L.faces);
vertices_R = double(struc_R.vertices);
faces_R = double(struc_R.faces);

subplot(1, 2, 1);
patch('Vertices', vertices_L, ...
      'Faces', faces_L, ...
      'FaceVertexCData', left.cdata, ...
      'FaceColor', 'interp', ...
      'EdgeColor', 'none');
view(2); axis equal off;
clim(limits);

subplot(1, 2, 2);
patch('Vertices', vertices_R, ...
      'Faces', faces_R, ...
      'FaceVertexCData', right.cdata, ...
      'FaceColor', 'interp', ...
      'EdgeColor', 'none');
view(2); axis equal off;
clim(limits);



end

