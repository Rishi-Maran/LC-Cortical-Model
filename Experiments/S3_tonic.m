%% Plot evoked response for different tonic strengths

clear; clc;
loadparam;

connparam.M = int32(0);

timeparam.T = timeparam.dt*ceil(0.6/timeparam.dt);

tau0 = coupparam.tau(1,4);

modparam.g(1) = 0.01;

tonic_array = [15, 20];

modparam.nu = 0; % tonic


title_array = cell(1, length(tonic_array));
for i = 1:length(tonic_array)
    title_array{i} = append('$Q_l = ', num2str(tonic_array(i)), '$');
end

fig = figure; fig.Position = [100 0 650 275]; hold;
t = tiledlayout(1, length(tonic_array), 'TileSpacing','loose');
% title(t, 'Effect of LC modulation on geometric propagation', 'Interpreter', 'latex', 'FontSize', 15);

% h_shift = 100;

for iter = 1:length(tonic_array)

    modparam.Qtonic = tonic_array(iter);
    ts = run_model_mex(timeparam,localparam,coupparam,connparam,modparam,stim,N_pop,N_nodes);
   
    %
    ax = nexttile(iter); hold on;
    set(ax, 'Layer', 'top');

    X = ts(1:N_nodes, :);
    X_max = max(X, [], 'all');
    X_norm = X / X_max;

    
    % Plot heatmap with time going to right
    imagesc(ax, 1000*(timeparam.dt*(1:ceil(0.6/timeparam.dt)) - stim.stimt(1)), 1:N_nodes, X_norm);  % transpose to make time vertical
    axis xy;      % ensures time goes down instead of up
    % set(ax, 'YTick', []);
    if iter == 1
        ylabel('Node', 'FontSize', 10, 'Interpreter', 'latex');
    end
    set(ax, 'TickLabelInterpreter', 'latex');
    xlim([0, 500]);
    ylim([1, N_nodes]);

    set(ax, 'XTick', 0:100:500);
    xlabel('$t$ (ms)', 'Interpreter', 'latex', 'FontSize', 10, 'Rotation', 0);
    
    title(title_array(iter), 'Interpreter', 'latex', 'FontSize', 10); 
    % yticks(0:50:N_nodes)
    yticks([]);
    colormap([linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)']); 
    
    % if iter == 2
    %     [~, idx] = max(ts(2*N_nodes+1, :));
    %     yl = yline(timeparam.dt*idx - stim.stimt(1) - tau0, 'Color', 'k', 'LineStyle', '--');
    %     if iter == 3
    %         yl.Label = "Time of Maximum Modulation";
    %     end
    % end

    % Plot stimuli position
    scatter(zeros(1, stim.stimnum), stim.pos(1:stim.stimnum), 5, [1, 0, 0], 'filled');


    clim([0, 1]);

    text(0.99, 0.99, sprintf('$\\phi_{\\mathrm{max}} = %.2f$',X_max), ...
        'Interpreter', 'latex', ...
        'Units','normalized', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','top', ...
        'FontSize',10);

    ax.Box = 'on';

    if iter == length(tonic_array)
        cb = colorbar(ax, 'southoutside');
        cb.TickLabelInterpreter = 'latex';
        cb.Layout.Tile = 'south';
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


end
%%


%
% save image
% saveas(gca, '2_geometric.svg');
exportgraphics(gcf, 'S2_tonic.tiff', 'Resolution', 600);
close(fig);