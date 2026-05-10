%% Effect of LC on wave propagation

clear; clc;
loadparam;

timeparam.T = 1;

modparam.g(1) = 0.01;

% fig = figure; fig.Position = [100 0 650 400]; hold;
% t = tiledlayout(1, 1, 'TileSpacing','loose');

% Timing
tau_c0 = modparam.tau_c;
mod_delay_array = 0:0.002:0.4;


max_activation_array = zeros(N_nodes, length(mod_delay_array));


parfor iter = 1:length(mod_delay_array)

    % Set nu_el (g(1))
    modparam1 = modparam;
    modparam1.tau_c = tau_c0 + mod_delay_array(iter);
    ts = run_multinode_sigmoid_mex(timeparam,localparam,coupparam,connparam,modparam1,stim,N_pop,N_nodes);

    X = ts(1:N_nodes, :);
    X_max = max(X, [], 2);
    
    max_activation_array(:, iter) = X_max;

    disp([mod_delay_array(iter), max(X_max)]);
    
end

save('S3_delay.mat', 'mod_delay_array', 'max_activation_array', 'timeparam', 'modparam');

%%

clear; clc;
loadparam;

load('S3_delay.mat');

fig = figure; fig.Position = [100 100 650 350]; hold;
t = tiledlayout(1, 1, 'TileSpacing', 'loose');

% Plot dissimilarity curves as heatmap
ax = gca; hold on;
set(ax,'YDir','normal','XDir','normal');
set(ax, 'Layer', 'top');
imagesc(1000*(mod_delay_array), 1:N_nodes, max_activation_array);
xlim([0, Inf]);
ylim([0, N_nodes]);
set(ax, 'TickLabelInterpreter', 'latex');
xlabel('Modulation Delay (ms)', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('Node', 'Interpreter', 'latex', 'FontSize', 10);
yticks([]);
colormap([linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)']); 
%Colorbar
cb = colorbar(ax,'southoutside'); 
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 8;
cb.Label.String = 'Peak Activation';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 10;

ax.Box = 'on';

%%
% save image
exportgraphics(gcf, 'S3_modulationdelay.tiff', 'Resolution', 600);
close(fig);
