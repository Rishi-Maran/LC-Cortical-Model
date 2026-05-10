%% Compute C(t) - effect of LC on non-local dissimilarity - for fine grained array of g values

clear; clc;
loadparam;

% Set stimulus position
stim.pos(1) = 51;

tau0 = coupparam.tau(1, 4);

g_array = [0:0.0002:0.05];

% Increase simulation time to 3;
timeparam.T = 3;
Nt = timeparam.T / timeparam.dt;

dissim_time_array = zeros(length(g_array), Nt);
dissim_bold_array = zeros(length(g_array), 1);

parfor iter = 1:length(g_array)

    modparam1 = modparam; modparam1.g(1) = g_array(iter); 

    ts0 = run_model_mex(timeparam,localparam,coupparam,connparam,modparam1,stim,N_pop,N_nodes);
    ts = run_model_mex(timeparam,localparam,coupparam,connparam_het,modparam1,stim,N_pop,N_nodes);

    dissim_time_array(iter, :) = compute_dissim(ts0(1:N_nodes, :), ts(1:N_nodes, :));

    bold0 = timeparam.dt*sum(ts0(1:N_nodes, :), 2)'; bold = timeparam.dt*sum(ts(1:N_nodes, :), 2)';

    % Check if either vector is all zeros
    if all(bold0 == 0) || all(bold == 0)
        dissim_bold_array(iter) = NaN;
    else
        dissim_bold_array(iter) = pdist2(bold0, bold, 'cosine');
    end

    disp([iter, max(dissim_time_array(iter, :)), dissim_bold_array(iter)]);

end
% save C(t) values in mat file for future plotting
save('5_dissim.mat', 'dissim_time_array', 'dissim_bold_array', 'timeparam', 'g_array');

%% Plot heatmap of C(t) and plot C_z for fine grained g value array

clear; clc;
loadparam;

load('5_dissim.mat');

num_subcols = 3;
num_cols = 3*num_subcols + 1;

fig = figure; fig.Position = [100 100 650 350]; hold;
t = tiledlayout(1, num_cols, 'TileSpacing', 'loose');

% Plot dissimilarity curves as heatmap
ax = nexttile(1, [1 num_subcols*2]); hold on;
set(ax,'YDir','normal','XDir','normal');
set(ax, 'Layer', 'top');
Nt = timeparam.T / timeparam.dt;
tau0 = coupparam.tau(1, 4);
imagesc(timeparam.dt*(1:Nt) - stim.stimt(1), g_array, dissim_time_array);
xlim([0, 1.2]);
ylim([min(g_array), max(g_array)]);
set(ax, 'TickLabelInterpreter', 'latex');
xlabel('$t$ (s)', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('$g$ ', 'Interpreter', 'latex', 'FontSize', 10, 'Rotation', 0);
colormap(flipud(pink)); 

%Colorbar
cb = colorbar(ax,'southoutside'); 
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 8;
cb.Label.String = '$C(t)$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 10;

ax.Box = 'on';

hold off;

ax = nexttile(2*num_subcols + 2, [1, num_subcols]); hold on;
set(ax,'YDir','normal','XDir','normal');
Nt = timeparam.T / timeparam.dt;
tau0 = coupparam.tau(1, 4);
plot(dissim_bold_array, g_array, 'LineWidth', 1, 'Color', 'k');
% xlim([0, 1.5]);
ylim([min(g_array), max(g_array)]);
set(ax, 'TickLabelInterpreter', 'latex');
xlabel('$C_z$', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('$g$ ', 'Interpreter', 'latex', 'FontSize', 10, 'Rotation', 0);

ax.Box = 'on';

hold off;
% Optional: Save as TIFF
% Uncomment the following lines if you want to save as TIFF
% exportgraphics(gcf, '5_dissim.tiff', 'Resolution', 600);
% close(fig);