%% Effect of gains on patchy cortical regions

%% Wave steering + Preferred stimuli

clear; clc;
loadparam;

timeparam.T = 0.75;

tau0 = coupparam.tau(1, 4);

% Set modulation gain
g_array = [0.02, 0.02, 0.02];

% Set stimuli
stimpos_array = zeros(3, 1);
stimpos_array(1) = 51;
stimpos_array(2) = 101;
stimpos_array(3) = 168;

N_array = [1; 2; 2];
ab_array = zeros(max(N_array), 2, 3);
ab_array(1:N_array(1), :, 1) = [51, 151];
ab_array(1:N_array(2), :, 2) = [101 168; 101 34];
ab_array(1:N_array(3), :, 3) = [101 168; 101 34];

% Set modulation region
kernel_array = zeros(3, N_nodes);
kernel_array(1, 49:53) = 1; kernel_array(1, 149:153) = 1;
kernel_array(2, 99:103) = 1; kernel_array(2, 166:170) = 1;
kernel_array(3, 99:103) = 1; kernel_array(3, 32:36) = 1;

num_subcols = 3;
num_cols = 3*num_subcols + 1;

fig = figure; fig.Position = [100 0 450 765]; hold;
t = tiledlayout(length(g_array), num_cols, 'TileSpacing','loose');

for iter = 1:3

    stim.pos(1) = stimpos_array(iter);

    modparam.g(1) = g_array(iter);
    modparam.kernel = kernel_array(iter, :);

    connparam.M = int32(2*N_array(iter));
    connparam.Chet(1:2*N_array(iter)) = 1e-3;
    connparam.ab(1:2:(2*N_array(iter)-1), :) = ab_array(1:N_array(iter), [1 2], iter);
    connparam.ab(2:2:(2*N_array(iter)), :) = ab_array(1:N_array(iter), [2 1], iter);

    ts = run_model_mex(timeparam,localparam,coupparam,connparam,modparam,stim,N_pop,N_nodes);
    %
    ax = nexttile(num_cols*(iter-1)+1, [1 num_subcols*2]); hold on;

    X = ts(1:N_nodes, :);
    X_max = max(X, [], 'all');
    X_norm = X / X_max;
    
    % Plot heatmap with time going to right
    imagesc((timeparam.dt:timeparam.dt:timeparam.T) - stim.stimt(1), 1:N_nodes, X_norm);  % transpose to make time vertical
    axis xy;      % ensures time goes down instead of up
    % set(ax, 'YTick', []);
    ylabel('Node', 'FontSize', 10, 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');
    xlim([0, 0.5]);
    ylim([1, N_nodes]);
    xlabel('$t$ (s)', 'Interpreter', 'latex', 'FontSize', 10, 'Rotation', 0);
    xticks(0:0.1:0.5);
    title('Hi', 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'none'); 
    yticks([]);
    colormap([linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)']); 

    % Plot modulation region
    idx_mod = find(kernel_array(iter, :) > 0);

    % Plot dummy scatter for legend
    plot(NaN, NaN, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '-', "LineWidth", 2);
    scatter(0.002*ones(size(idx_mod)), idx_mod, 10, [0.5, 0.5, 0.5], 'filled', 'square', 'MarkerFaceAlpha', 1);
    % Plot stimuli position
    scatter(zeros(1, stim.stimnum), stim.pos(1:stim.stimnum), 5, [1, 0, 0], 'filled', 'MarkerFaceAlpha', 1);

    clim([0, 1]);

    text(0.99, 0.99, sprintf('$\\phi_{\\mathrm{max}} = %.2f$',X_max), ...
        'Interpreter', 'latex', ...
        'Units','normalized', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','top', ...
        'FontSize',10);
    
    ax.Box = 'on';

    cb = colorbar(ax, 'southoutside');
    cb.TickLabelInterpreter = 'latex';
    % cb.Layout.Tile = 'south';
    cb.Ticks = [0 1];
    cb.TickLabels = {'$0$', '$\phi_{\mathrm{max}}$'}; 
    cb.FontSize = 8;
    cb.Label.String = '$\phi(t)$';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = 10;
    
    % Move label up
    pos = cb.Label.Position;
    pos(2) = pos(2) + 0.8;   % adjust y offset
    cb.Label.Position = pos;

    % Draw lrc positions
    yline(reshape(ab_array(:, :, iter), 1, []), 'Color', [0, 0.5, 0], 'LineStyle', '--', 'LineWidth', 1);

    hold off;
    
    pause(0.5);

    ax = nexttile(num_cols*(iter - 1) + 2*num_subcols + 2, [1, num_subcols]); hold on;

    axis xy;      % ensures time goes down instead of up
    xlim([0, 1]);
    xlabel('Peak Activation', 'Interpreter', 'latex', 'FontSize', 10);
    ax.XTick = [0 1];
    ax.XTickLabel = {'0' '$\phi_{\mathrm{max}}$'};

    set(ax, 'YTick', []);
    ylim([1, N_nodes]);
    set(ax, 'TickLabelInterpreter', 'latex');
    title('Hi', 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'none'); 
    % Plot heatmap with time going down
    plot(max(X_norm, [], 2), 1:N_nodes, 'Color', 'r', 'LineWidth', 1);  % transpose to make time vertical
    
    set(ax, 'YTick', []);
    % Plot stimuli position
    scatter(zeros(1, stim.stimnum), stim.pos(1:stim.stimnum), 5, [1, 0, 0], 'filled');

    pause(0.5);

end 

for iter = 1
    ax = nexttile(num_cols*(iter-1)+1, [1 num_subcols*2]); hold on;
    h = legend(ax, {"Modulated Region", "", "", "", "", "", "", "", ""}, 'Location', 'east', 'Interpreter', 'latex', 'FontSize', 8, 'Box', 'off');  
    hold off;
end

%

% Optional: Save as TIFF
% Uncomment the following lines if you want to save as TIFF
% exportgraphics(gcf, '6_patchy_lrc.tiff', 'Resolution', 600);
% close(fig);
