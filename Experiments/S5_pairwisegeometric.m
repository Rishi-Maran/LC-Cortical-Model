%% Max distance versus g - endpoint modulation

clear; clc;
loadparam;

% Change stimulus position
stim.pos(1) = lrc_a;

timeparam.T = 0.5;

N = double(N_nodes);

min_g = 0; min_dist = 1;
max_g = 1; max_dist = 100;
d_g = 0.001; d_dist = 1;

threshold_array = [5, 10, 50, 100];

contour_points = cell(length(threshold_array), 1);

parfor iter = 1:length(threshold_array)

    threshold = threshold_array(iter);

    base = zeros(1, N);

    modparam1 = modparam;

    g = min_g; dist = min_dist;
    points = [];
    fval = f_sim(g, dist, modparam1, lrc_a, base, timeparam, localparam, ...
                 coupparam, connparam, stim, N_pop, N_nodes);

    disp([threshold, g, dist, fval]);
    
    while g < max_g && dist < max_dist
    
        % --- march along g at current dist ---
        g_next = g;
        while g_next < max_g
            g_next = g_next + d_g;
            f_next = f_sim(g_next, dist, modparam1, lrc_a, base, timeparam, localparam, ...
                 coupparam, connparam, stim, N_pop, N_nodes);
            disp([threshold, g_next, dist, f_next]);
            if (fval - threshold)*(f_next - threshold) < 0
                % Linear interpolation for threshold crossing
                g_cross = g + (threshold - fval)*(g_next - g)/(f_next - fval);
                points = [points; g_cross, dist];
                break  % stop marching in g, now step dist
            end
        end

        if g_next == max_g
            break;
        end

        g = g_next;
        fval = f_next;
    
        % --- march along dist at current g ---
        dist_next = dist;
        while dist_next < max_dist
            dist_next = dist_next + d_dist;
            f_next = f_sim(g, dist_next, modparam1, lrc_a, base, timeparam, localparam, ...
                 coupparam, connparam, stim, N_pop, N_nodes);
            disp([threshold, g, dist_next, f_next]);
            if (fval - threshold)*(f_next - threshold) < 0
                % Linear interpolation for threshold crossing
                dist_cross = dist + (threshold - fval)*(dist_next - dist)/(f_next - fval);
                points = [points; g, dist_cross];
                break  % stop marching in dist, now step g
            end
            dist = dist_next;
            fval = f_next;
        end

        if dist_next == max_dist
            break;
        end

        dist = dist_next;
        fval = f_next;
    
    end

    contour_points{iter} = points;

end


% Save results
save('contour_points.mat', 'contour_points', 'threshold_array', 'min_g', 'max_g', 'd_g', 'min_dist', 'max_dist', 'd_dist');

% --- Helper function (subfunction, outside script for parfor) ---
% function fval = f_sim(g, dist, modparam1, lrc_a, base, timeparam, localparam, ...
%                        coupparam, connparam, stim, N_pop, N_nodes)
%     modparam1.g(1) = g;
%     lrc_b = lrc_a + dist;
%     base(:) = 0;
%     base([lrc_a, lrc_b]) = 1;
%     modparam1.kernel = base;
% 
%     ts = run_model_mex(timeparam, localparam, coupparam, connparam, ...
%                                    modparam1, stim, N_pop, N_nodes);
%     fval = max(ts(lrc_b, :));
% end

%% Plot contour points

clear; clc;
loadparam;

load('contour_points.mat')


fig = figure; fig.Position = [100, 0, 650, 300]; hold;
l = cell(4, 1);
for iter = 1:4
    hi = contour_points{iter};
    plot(hi(:,2), hi(:,1), 'LineWidth', 1);
    l{iter} = append('Target Activation $\approx$ ', num2str(threshold_array(iter)));
end

legend(l, 'Interpreter', 'latex', 'Location', 'northwest');
xlabel('Stimulus-Target Distance (D)', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$g$', 'Interpreter', 'latex', 'Rotation', 0, 'FontSize', 15);
ax = gca; ax.TickLabelInterpreter = 'latex';

% Optional: Save as TIFF
% Uncomment the following lines if you want to save as TIFF
% exportgraphics(gcf, 'S4_contour.tiff', 'Resolution', 600);
% close(fig);
