%% Effect of LC timing on wave propagation

clear; clc;
loadparam;

connparam.M = int32(0);

timeparam.T = 0.75;

modparam.g(1) = 0.01;

fig = figure; fig.Position = [100 0 650 495]; hold;
t = tiledlayout(2, 2, 'TileSpacing','loose');

% Timing specificity
tau_c0 = modparam.tau_c;
mod_delay_array = [0.025, 0.05];

title_array = cell(1, 2);
title_array{1} = append('$', num2str(1000*mod_delay_array(1)), '\mathrm{ms}$ modulation delay');
title_array{2} = append('$', num2str(1000*mod_delay_array(2)), '\mathrm{ms}$ modulation delay');

for iter = 1:length(mod_delay_array)

    % Set timing of LC firing
    modparam.tau_c = tau_c0 + mod_delay_array(iter);

    ts = run_model_mex(timeparam,localparam,coupparam,connparam,modparam,stim,N_pop,N_nodes);
    Nt_new = ceil((0.5 + stim.stimt(1)) / timeparam.dt);
    ts = ts(:, 1:Nt_new);
    X = ts(1:N_nodes, :);
    X_max = max(X, [], 'all');
    X_norm = X / X_max;
    
    ax0(iter) = nexttile(t); hold on;
    set(ax0(iter), 'Layer', 'top');
    hold(ax0(iter), 'on');
   
    % Plot heatmap with time going to right
    imagesc(ax0(iter), 1000*(timeparam.dt*(1:Nt_new) - stim.stimt(1)), 1:N_nodes, X_norm);  % transpose to make time vertical
    axis(ax0(iter), 'xy');      % ensures time goes down instead of up
    ylabel(ax0(iter), 'Node', 'FontSize', 10, 'Interpreter', 'latex');
    set(ax0(iter), 'TickLabelInterpreter', 'latex');
    xlim(ax0(iter), [0, 500]);
    ylim(ax0(iter), [1, N_nodes]);

    set(ax0(iter), 'XTick', 0:100:500);  % Always show tick marks
    set(ax0(iter), 'XTickLabel', []);  % Hide labels except on last row

    title(title_array(iter), 'Interpreter', 'latex', 'FontSize', 10); 
    % yticks(0:50:N_nodes)
    yticks(ax0(iter), []);
    colormap(ax0(iter), [linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)']); 
    ax0(iter).Box = 'on';

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
    plot(ax0(iter), 1000*(timeparam.dt*(1:Nt_new) - stim.stimt(1)), 1 + ts(2*N_nodes+1, :), 'Color', [0.5, 0.5, 0.5]);
    yline(1, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1, 'LineStyle', '--', 'Alpha', 0.5);
    ax0(iter).YAxis(2).Color = [0.5 0.5 0.5];   % right y-axis color
    ylim([0, 12]);
    if iter == 1
        yticks([0, 1, 2, 3]);
        ax0(iter).YAxis(2).FontSize = 8;
    else
        yticks([]);
    end

    text(0.99, 1/12, 'Modulation Factor', ...
        'Interpreter', 'latex', ...
        'Units','normalized', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','bottom', ...
        'Color', [0.5 0.5 0.5], ...
        'FontSize',8); 
 
    ax0(iter).Box = 'on';


    xl = xline(1000*coupparam.tau(1, 4), 'Interpreter', 'latex', ...
        'FontSize', 8, 'Color', 'k', 'LabelOrientation', 'horizontal', ...
        'LineStyle', '--');
    if iter == 1
        xl.Label = ['Time of', newline, 'Activation'];
    end

    hold off;

    pause(0.5);

end

% Spatial specificity

% Reset parameters
loadparam;
connparam.M = int32(0);
timeparam.T = 0.75;
modparam.g(1) = 0.01;

N = double(N_nodes);

num_trials = 2;

% Modulation kernel generation
kernel_array = zeros(num_trials, N_nodes);
kernel_array(1, 1 + 0:round(6*N/10)) = 1;
kernel_array(2, 1 + [0:round(3*N/10) round(4*N/10):round(6*N/10)]) = 1;

for iter = 1:num_trials

    % Set cortical targets of LC
    modparam.kernel = kernel_array(iter, :);

    ts = run_model(timeparam,localparam,coupparam,connparam,modparam,stim,N_pop,N_nodes);
    X = ts(1:N_nodes, :);
    X_max = max(X, [], 'all');
    X_norm = X / X_max;

    ax1(iter) = nexttile(t);  % or ax1(iter) = axes(t); if you're indexing manually
    set(ax1(iter), 'Layer', 'top');
    hold(ax1(iter), 'on');
    
    % Plot your main heatmap
    imagesc(ax1(iter), 1000*((timeparam.dt:timeparam.dt:timeparam.T) - stim.stimt(1)), 1:N_nodes, X_norm);
    axis(ax1(iter), 'xy');
    ylabel(ax1(iter), 'Node', 'FontSize', 10, 'Interpreter', 'latex');
    set(ax1(iter), 'TickLabelInterpreter', 'latex');
    xlim(ax1(iter), [0, 500]);
    ylim(ax1(iter), [1, N_nodes]);

    set(ax1(iter), 'XTick', 0:100:500);  % Always show tick marks
    xlabel('$t$ (ms)', 'Interpreter', 'latex', 'FontSize', 10, 'Rotation', 0);

    yticks(ax1(iter), []);
    colormap(ax1(iter), [linspace(1, 1, 256)' linspace(1, 0, 256)' linspace(1, 0, 256)']);
    ax1(iter).Box = 'on';

    % Plot target nodes
    k_indices = find(kernel_array(iter, :) > 0);  % get indices where arr == 1
    x = zeros(size(k_indices));  % x = 0 for all points
    y = k_indices;
    scatter(ax1(iter), x, y, 10, [0.5, 0.5, 0.5], 'filled');  % plot filled circles

    % Plot stimuli position
    scatter(ax1(iter), zeros(1, stim.stimnum), stim.pos(1:stim.stimnum), 5, [1, 0, 0], 'filled');

    text(0.02, 0.02, 'Modulated Region', ...
        'Interpreter', 'latex', ...
        'Units', 'normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'Color', [0.5 0.5 0.5], ...
        'FontSize',8); 

    clim(ax1(iter), [0 1]);

    text(0.99, 0.99, sprintf('$\\phi_{\\mathrm{max}} = %.2f$',X_max), ...
        'Interpreter', 'latex', ...
        'Units','normalized', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','top', ...
        'FontSize',10);


    if iter == num_trials
        cb = colorbar(ax1(iter), 'southoutside');
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
    hold(ax1(iter), 'off');

    pause(0.5);

end

drawnow;
%
% Optional: Save as TIFF
% Uncomment the following lines if you want to save as TIFF
% exportgraphics(gcf, '4_variability.tiff', 'Resolution', 600);
% close(fig);