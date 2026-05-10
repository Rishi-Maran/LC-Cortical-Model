%% Effect of LC with varying modulation strengths

%% No LRCs

clear; clc;
loadparam_dti;

connparam.M = int32(0);

timeparam.T = 0.5;

tau0 = coupparam.tau(1, 4);

% g_array = 0:0.0004:0.0008;
g_array = [0.00, 0.01];

fig = figure; fig.Position = [100 100 650 500]; hold;
t = tiledlayout(3, 3*length(g_array));
% title(t, 'Effect of LC Modulation on Evoked Response', 'Interpreter', 'latex', 'FontSize', 15);

% Set coupling from stimulus to LC
modparam.nu = coupparam.nu(4, 5);

for iter = 1:length(g_array)

    modparam.g = g_array(iter)*[1; 0; 0; 0; 0];
    ts = run_model_mex(timeparam,localparam,coupparam,connparam,modparam,stim,N_pop,N_nodes);
    ts = ts(1:N_nodes, :);
    ax = nexttile(t, 1 + (iter-1)*3, [3, 2]); hold on;

    for i = 1:N_nodes
        [~, idx] = max(ts(i, :));
        peak_time_array(i) = timeparam.dt*idx - stim.stimt(1) - tau0;
    end
s
    [~, sort_idx] = sort(peak_time_array, 'ascend');

    % X = ts(sort_idx, :) ./ max(ts(sort_idx, :), [], 2);
    X = ts(sort_idx, :);
    X_max = max(X, [], 'all');
    X_norm = X / X_max;

    set(ax, 'Layer', 'top');
    % Plot heatmap with time going down
    imagesc(1:N_nodes, 1000*((timeparam.dt:timeparam.dt:timeparam.T) - stim.stimt(1)), X_norm');  % transpose to make time vertical
    axis ij;      % ensures time goes down
    xlim([0.5, N_nodes + 0.5]);
    ylim([0 280]);
    set(ax, 'XTick', []);
    ax.XAxis.MinorTick = 'on';
    set(ax, 'TickLabelInterpreter', 'latex');

    if iter == 1
        ylabel('$t$ (ms)', 'Interpreter', 'latex', 'FontSize', 10, 'Rotation', 90);
        title('No Modulation ($g = 0$)', 'Interpreter', 'latex', 'FontSize', 10); 
        ax.XAxis.MinorTickValues = 0:100:300;
    else
        set(ax, 'YTickLabel', []);
        title(append('Phasic Modulation ($g = ',num2str(g_array(iter)), '$)'), 'Interpreter', 'latex', 'FontSize', 10);
    end
    xlabel('Node', 'FontSize', 10, 'Interpreter', 'latex');
    
    colormap([linspace(1, 1, 256)' linspace(1, 0, 256)', linspace(1, 0, 256)']); 
    % xlim([0.5, N_nodes - 0.5]);
    % clim(max(X, [], 'all')*[0, 1]);
    % yl = yline(tau0, 'Color', [0.5, 0.5, 0.5]);
    % if iter == 1
    %     yl.Label = "Time of Cortex Activation";
    %     yl.Interpreter = "latex"; yl.FontSize = 8;
    %     yl.LabelVerticalAlignment = "bottom";
    % end

    % scatter(1, 0, 4, 'r', 'filled');

    text(0.99, 0.99, sprintf('$\\phi_{\\mathrm{max}} = %.2f$',X_max), ...
    'Interpreter', 'latex', ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top', ...
    'FontSize',10);


    cb = colorbar(ax, 'southoutside', 'TickLabelInterpreter', 'latex');
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

    % Plot with incresing time to peak
    % plot(1:N_nodes, peak_time_array(sort_idx), 'LineWidth', 0.5, 'Color', 'k');

    % if iter == 1
    %     yPos = 0.35;
    %     xPos = 110;
    %     text(ax, xPos, yPos, sprintf('Time of\nPeak Activation'), ...
    %         'Interpreter', 'latex', ...
    %         'HorizontalAlignment','center', ...
    %         'VerticalAlignment','middle', ...
    %         'FontSize',8,'FontWeight','bold');
    % end

    % pos = ax.Position;
    % annotation('rectangle', pos, 'Color', 'black', 'LineWidth', 0.01);
    ax.Box = 'on';

    hold off;

    times = [0.05, 0.1, 0.2];
 
    for i = 1:3
        ax1 = nexttile(t, 3 + (iter-1)*3 + 3*(i-1)*length(g_array)); hold on;
        axis off; title(sprintf('$t = %.0f$ ms', 1000*times(i)), 'Interpreter', 'latex', 'FontSize', 8);
        idx = round((stim.stimt(1) + times(i))/timeparam.dt);
        surf_parcel_snapshot_midthickness_left(ts(:, idx)/X_max,[],[0, 1],[]);
    end
    pause(0.5);

    % nexttile(t, 3 + (iter-1)*3); axis off; title('$t = 0.05$', 'Interpreter', 'latex', 'FontSize', 8);
    % nexttile(t, 3 + (iter-1)*3 + 3*length(g_array)); axis off; title('$t = 0.1$', 'Interpreter', 'latex', 'FontSize', 8);
    % nexttile(t, 3 + (iter-1)*3 + 6*length(g_array)); axis off; title('$t = 0.2$', 'Interpreter', 'latex', 'FontSize', 8);

    pause(0.5);


end

%
% save image
exportgraphics(gcf, 'S1_response_dti.tiff', 'Resolution', 300);
close(fig);

%% Plot 3D snapshots of evoked response by LC modulation strength

clear; clc;
loadparam_dti;

connparam.M = int32(0);

% timeparam.T = 1;
timeparam.T = round(0.35/timeparam.dt)*timeparam.dt;

g_array = 0:0.0025:0.005;

% Set coupling from stimulus to LC
modparam.nu = coupparam.nu(4, 5);

for iter = 1:length(g_array)

    modparam.g = g_array(iter)*[1; 0; 0; 0; 0];
    ts = run_multinode_sigmoid(timeparam,localparam,coupparam,connparam,modparam,stim,N_pop,N_nodes);
    ts = ts(1:N_nodes, :);

    times = [0.05, 0.1, 0.2];
   
    for i = 1:3
        idx = round((stim.stimt(1) + times(i))/timeparam.dt);
        surf_parcel_snapshot_midthickness_left(ts(:, idx),stim.pos(1),max(ts, [], 'all')*[0, 1],[]);
    end
    pause(0.5);


end

%% Simulate evoked response with and without gain

clear; clc;
loadparam_dti;

connparam.M = int32(0);

% timeparam.T = 1;
timeparam.T = round(0.35/timeparam.dt)*timeparam.dt;

g_array = 0:0.005:0.01;

% Set coupling from stimulus to LC
modparam.nu = coupparam.nu(4, 5);

for iter = 1:length(g_array)

    modparam.g = g_array(iter)*[1; 0; 0; 0; 0];
    ts = run_multinode_sigmoid(timeparam,localparam,coupparam,connparam,modparam,stim,N_pop,N_nodes);
    ts = ts(1:N_nodes, round(stim.stimt(1)/timeparam.dt):end);

    surf_parcel_animate_midthickness_left(ts,timeparam.dt, ...
        0,timeparam.T,100,stim.pos(1),[],sprintf('4_no_lrcs_g%.1e.gif', g_array(iter)),...
        sprintf('g = %.1e, No LRCs', g_array(iter)),[]);

    pause(0.5);


end

%% Simulate evoked response with LRCs, with and without gain

clear; clc;
loadparam_dti;

timeparam.T = round(0.35/timeparam.dt)*timeparam.dt;

g_array = 0:0.005:0.01;

% Set coupling from stimulus to LC
modparam.nu = coupparam.nu(4, 5);

for iter = 1:length(g_array)

    modparam.g = g_array(iter)*[1; 0; 0; 0; 0];
    ts = run_multinode_sigmoid(timeparam,localparam,coupparam,connparam,modparam,stim,N_pop,N_nodes);
    ts = ts(1:N_nodes, round(stim.stimt(1)/timeparam.dt):end);

    surf_parcel_animate_midthickness_left(ts,timeparam.dt, ...
        0,Inf,100,stim.pos(1),[],sprintf('4_g%.1e.gif', g_array(iter)), ...
        sprintf('g = %.1e, Infinite speed LRCs', g_array(iter)),connparam.ab(1:N_lrcs, :));

    pause(0.5);


end

%% Simulate evoked response at most posterior point with LRCs, with and without gain

clear; clc;
loadparam_dti;

% Determine parcel of most posterior vertex of mesh
v_data = gifti('../Model/Data/fs_LR.32k.L.midthickness.surf.gii').vertices;
[~, v_idx] = min(v_data(:, 2)); clear v_data;
V_P_map = dlmread('../Model/Data/fsLR_32k_Schaefer400-lh.txt');
stim.pos(1) = V_P_map(v_idx); clear V_P_map;

timeparam.T = round(0.35/timeparam.dt)*timeparam.dt;

g_array = 0:0.01:0.04;

% Set coupling from stimulus to LC
modparam.nu = coupparam.nu(4, 5);

for iter = 1:length(g_array)

    modparam.g = g_array(iter)*[1; -1; 0; 0; 0];
    ts = run_multinode_sigmoid(timeparam,localparam,coupparam,connparam,modparam,stim,N_pop,N_nodes);
    ts = ts(1:N_nodes, round(stim.stimt(1)/timeparam.dt):end);

    surf_parcel_animate_midthickness_left(ts,timeparam.dt, ...
        0,Inf,100,stim.pos(1),[],0, ...%sprintf('4_g%.2f_posterior.gif', g_array(iter)), ...
        sprintf('g = %.2f, Infinite speed LRCs', g_array(iter)),connparam.ab(1:N_lrcs, :));

    pause(0.5);


end


