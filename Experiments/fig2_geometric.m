%% Effect of LC with varying modulation strengths

%% Effect of LC on wave propagation

clear; clc;
loadparam;

connparam.M = int32(0);

timeparam.T = 0.75;

tau0 = coupparam.tau(1,4);

g_array = [0, 0.01, 0.01];
tonic_array = zeros(size(g_array)); tonic_array(3) = 10;
mod_nu_array = modparam.nu*ones(size(g_array)); mod_nu_array(3) = 0;

title_array = cell(1, 3);
title_array{1} = append('No Modulation ($g = 0$)');
title_array{2} = 'Phasic Modulation';
title_array{3} = 'Tonic Modulation';

num_subcols = 3;
num_cols = 3*num_subcols + 1;

fig = figure; fig.Position = [100 0 450 650]; hold;
t = tiledlayout(length(g_array), num_cols, 'TileSpacing','loose');
% title(t, 'Effect of LC modulation on geometric propagation', 'Interpreter', 'latex', 'FontSize', 15);

% h_shift = 100;

for iter = 1:length(g_array)

    % Set nu_el (g(1))
    modparam.g(1) = g_array(iter); 
    modparam.Qtonic = tonic_array(iter);
    modparam.nu = mod_nu_array(iter);
    ts = run_model_mex(timeparam,localparam,coupparam,connparam,modparam,stim,N_pop,N_nodes);
    %
    ax = nexttile(num_cols*(iter-1)+1, [1 num_subcols*2]); hold on;
    set(ax, 'Layer', 'top');

    X = ts(1:N_nodes, :);
    X_max = max(X, [], 'all');
    X_norm = X / X_max;

    
    % Plot heatmap with time going to right
    imagesc(ax, 1000*((timeparam.dt:timeparam.dt:timeparam.T) - stim.stimt(1)), 1:N_nodes, X_norm);  % transpose to make time vertical
    axis xy;      % ensures time goes down instead of up
    % set(ax, 'YTick', []);
    ylabel('Node', 'FontSize', 10, 'Interpreter', 'latex');
    set(ax, 'TickLabelInterpreter', 'latex');
    xlim([0, 500]);
    ylim([1, N_nodes]);

    set(ax, 'XTick', 0:100:500);
    if iter < length(g_array)
        set(ax, 'XTickLabel', []);
    else
        xlabel('$t$ (ms)', 'Interpreter', 'latex', 'FontSize', 10, 'Rotation', 0);
    end
    title(title_array(iter), 'Interpreter', 'latex', 'FontSize', 10); 
    % yticks(0:50:N_nodes)
    yticks([]);
    colormap([linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)']); 

    % Plot stimuli position
    scatter(zeros(1, stim.stimnum), stim.pos(1:stim.stimnum), 5, [1, 0, 0], 'filled');


    clim([0, 1]);

    text(0.99, 0.99, sprintf('$\\phi_{\\mathrm{max}} = %.2f$',X_max), ...
        'Interpreter', 'latex', ...
        'Units','normalized', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','top', ...
        'FontSize',10);

    yyaxis right;
    plot(ax, 1000*((timeparam.dt:timeparam.dt:timeparam.T) - stim.stimt(1)), 1 + ts(2*N_nodes+1, :), 'Color', [0.5, 0.5, 0.5]);
    yline(1, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1, 'LineStyle', '--', 'Alpha', 0.5);
    % ylabel('Modulation Factor', 'Rotation', -90, 'Color', [0.5, 0.5, 0.5], 'Interpreter', 'latex');
    % ax = gca;  
    ax.YAxis(2).Color = [0.5 0.5 0.5];   % right y-axis color
    ylim([0, 12]);
    yticks([0, 1, 2, 3]);
    ax.YAxis(2).FontSize = 8;
    % ytickformat('percentage');
    if iter == 1 || iter == 3
        Text = append('Modulation Factor = $', num2str(1 + ts(2*N_nodes+1, end)), '$');
    else
        Text = 'Modulation Factor';
    end
    text(0.99, 1/12, Text, ...
        'Interpreter', 'latex', ...
        'Units','normalized', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom', ...
        'Color', [0.5 0.5 0.5], ...
        'FontSize',8);    
 
    ax.Box = 'on';

    if iter == length(g_array)
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
    end
    
    hold off;

    pause(0.5);

    ax = nexttile(num_cols*(iter - 1) + 2*num_subcols + 2, [1, num_subcols]); hold on;

    axis xy;      % ensures time goes down instead of up
    xlim([0, 1]);
    if iter == length(g_array)
        xlabel('Peak Activation', 'Interpreter', 'latex', 'FontSize', 10);
        ax.XTick = [0 1];
        ax.XTickLabel = {'0' '$\phi_{\mathrm{max}}$'};
    else
        ax.XTick = [];
    end

    set(ax, 'YTick', []);
    ylabel('Node', 'FontSize', 10, 'Interpreter', 'latex');
    ylim([1, N_nodes]);
    set(ax, 'TickLabelInterpreter', 'latex');
    title('Hi', 'Interpreter', 'latex', 'FontSize', 10, 'Color', 'none'); 
    % Plot heatmap with time going down
    plot(max(X_norm, [], 2), 1:N_nodes, 'Color', 'r', 'LineWidth', 1);  % transpose to make time vertical
    set(ax, 'YTick', []);
    pause(0.5);
end

% Optional: Save as TIFF
% Uncomment the following lines if you want to save as TIFF
% exportgraphics(gcf, '2_geometric.tiff', 'Resolution', 600);
% close(fig);