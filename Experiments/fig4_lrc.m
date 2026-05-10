%% Geometric propagation with increasing gain

clear; clc;
loadparam;

% set lrc position
a = 51; b = 151;

timeparam.T = 0.75;

connparam.M = int32(0);

tau0 = coupparam.tau(1,4);

% Change stimulus position
stim.pos(1) = 51;

g_array = [0.01, 0.02, 0.03];

title_array = cell(1, 3);
for iter = 1:3
    title_array{iter} = append('$g = ', num2str(g_array(iter)), '$');
end

fig = figure; fig.Position = [100 0 1000 495]; hold;
t = tiledlayout(2, length(g_array), 'TileSpacing','loose');

for iter1 = 1:2

    for iter = 1:length(g_array)

        modparam.g(1) = g_array(iter);

        if iter1 == 1
            ts = run_model_mex(timeparam,localparam,coupparam,connparam,modparam,stim,N_pop,N_nodes);
        else
            ts = run_model_mex(timeparam,localparam,coupparam,connparam_het,modparam,stim,N_pop,N_nodes);
        end
        %
        ax = nexttile; hold on;
        set(ax, 'Layer', 'top');
        
        X = ts(1:N_nodes, :);
        X_max = max(X, [], 'all');
        X_norm = X / X_max;
    
        
        % Plot heatmap with time going to right
        h = imagesc(1000*((timeparam.dt:timeparam.dt:timeparam.T) - stim.stimt(1)), 1:N_nodes, X_norm);  % transpose to make time vertical
        axis xy;      % ensures time goes down instead of up
        set(ax, 'TickLabelInterpreter', 'latex');
        xlim([0, 500]);
        ylim([1, N_nodes]);
        set(ax, 'XTick', 0:100:500);  % Always show tick marks
        
        if iter1 == 1
            set(ax, 'XTickLabel', []);  % Hide labels except on last row
        else
            xlabel('$t$ (ms)', 'Interpreter', 'latex', 'FontSize', 10, 'Rotation', 0);
        end

        ylabel('Node', 'Interpreter', 'latex', 'FontSize', 10);
        title(title_array(iter), 'Interpreter', 'latex', 'FontSize', 10);
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

        if iter1 == 2    
            % Draw dotted lines for LRC
            yline(a, 'Color', [0.5, 0.5, 1], 'LineStyle', '--', 'LineWidth', 0.5);
            yline(b, 'Color', [0.5, 0.5, 1], 'LineStyle', '--', 'LineWidth', 0.5);
        end
    
        ax.Box = 'on';
    
        if iter == length(g_array) && iter1 == 2
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

end

%

% Optional: Save as TIFF
% Uncomment the following lines if you want to save as TIFF
% exportgraphics(gcf, '4_lrc.tiff', 'Resolution', 600);
% close(fig);

%% Effect of LC on non-local dissimilarity - C(t)

clear; clc;
loadparam;

% set lrc position
a = 51; b = 151;

% Set stimulus position
stim.pos(1) = 51;

timeparam.T = 1.75;

connparam.M = int32(0);

tau0 = coupparam.tau(1, 4);

g_array = [0.01, 0.02, 0.03];

Nt = timeparam.T / timeparam.dt;

fig = figure; fig.Position = [100 100 1000 300]; hold;
xlim([0, 1500]);

for iter = 1:length(g_array)

    modparam.g(1) = g_array(iter); 

    % Simulate responses
    ts0 = run_model_mex(timeparam,localparam,coupparam,connparam,modparam,stim,N_pop,N_nodes);
    ts = run_model_mex(timeparam,localparam,coupparam,connparam_het,modparam,stim,N_pop,N_nodes);

    dissim_time = compute_dissim(ts0(1:N_nodes, :), ts(1:N_nodes, :));
    plot(1000*(timeparam.dt*(1:Nt) - stim.stimt(1)), dissim_time, 'LineWidth', 1);

    l{length(g_array) - iter + 1} = append('$g = ', num2str(g_array(iter)), '$');
    
    pause(0.5);

end
%
% Reorder layering of curves, put mean curves on top in with least N value curve on top
h = get(gca,'Children');
set(gca, 'Children', [h(length(g_array):-1:1)]);


xlabel("$t \ (\mathrm{ms})$", 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$C(t)$', 'Interpreter', 'latex', 'Rotation', 0, 'FontSize', 12)
%
l = legend(l, 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 10);

ax = gca;
set(ax, 'TickLabelInterpreter', 'latex');
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0:100:1500;   % R2015b+ only

%

% Optional: Save as TIFF
% Uncomment the following lines if you want to save as TIFF
% exportgraphics(gcf, '4_dissim_curve.tiff', 'Resolution', 1000);
% close(fig);


