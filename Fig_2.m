%% FIGURE 2 | MOVEMENT-BASED TRAP-CROPPING MODEL
% Paper: Insect movement and spatial dynamics determine trap crop effectiveness
% Journal: Ecological Modelling (2026)
%
% PURPOSE
%   Simulate one-dimensional insect movement in a cash-crop / trap-crop
%   landscape with temporary emigration, and reproduce the four-panel
%   presentation of Figure 2:
%     (a) Position distribution at day 0
%     (b) Position distribution at day 2
%     (c) Position distribution at day 30
%     (d) Percentages in four mutually exclusive population states
%
% MODEL OVERVIEW
%   Cash crop:         0 <= x <= k*L
%   Trap crop:       k*L <  x <= L
%   Emigration:        x < 0 or x > L (temporary; re-entry permitted)
%   Removed:           insect permanently removed from the simulation
%
%   During each movement step, insect i updates its position as
%       x_i(t+1) = x_i(t) + delta_i(t+1) * s_i(t+1),
%   where s_i is half-normal with scale sigma in the cash crop and outside
%   the field. Inside the trap crop, the scale is rho*sigma.
%
%   Cash-crop directions are chosen by a biased correlated random walk:
%     * With probability 1-omega, use directional persistence alpha.
%     * With probability omega, use spatial attraction beta(x).
%       beta(x) = beta_max * (x/(k*L))^nu, for 0 <= x <= k*L.
%   Trap-crop and emigrated insects independently choose direction +1 or
%   -1 with probability 1/2 at every movement step.
%
%   Removal mechanisms:
%     * On entering the trap crop from the cash crop, remove with
%       probability gamma_entry (each new entry is another opportunity).
%     * During residence inside the trap, remove at each time step with
%       probability gamma_dt = 1 - (1-gamma_day)^(dt/24), for dt in hours.
%
% REPRODUCIBILITY AND VALIDATION
%   MATLAB R2016a-compatible; no external toolboxes are required.
%   The random seed is fixed for repeatability on compatible MATLAB systems.
%   Four population categories are mutually exclusive and exhaustive at each
%   recorded time. The script checks this conservation identity explicitly.
%
% FILE OUTPUT
%   Automatically export a tightly cropped PNG in addition to displaying
%   the figure. Set save_figure to false to disable PNG export.

clear; clc; close all;
rng(1, 'twister');                % Fixed random seed for reproducibility

%% 1. BASELINE PARAMETERS
% Time and population
dt           = 0.5;              % Duration of one step (hours)
t_days       = 30;               % Total simulated time (days)
Npart        = 5000;             % Initial number of insects
p_steps      = round(t_days * 24 / dt) + 1; % Include t = 0

% Spatial layout and movement
L            = 100;              % Total field length (m)
k_cash       = 0.80;             % Cash-crop fraction, k in the manuscript
zone_left    = k_cash * L;       % Cash / trap interface at x = 80 m
sigma_cash   = 1;                % Half-normal step-length scale (m)
rho          = 0.1;              % Trap-step scale factor (0 < rho <= 1)
sigma_trap   = rho * sigma_cash; % Arrested movement scale (m)

% Cash-crop BCRW
alpha        = 0.2;              % Directional persistence (-1 <= alpha <= 1)
omega        = 0.5;              % Chance of selecting cue-driven movement
beta_max     = 0.1;              % Maximum directional bias at x = k*L
nu           = 1;                % Spatial attraction exponent

% Trap-crop removal
gamma_entry  = 0.1;              % Probability of removal at each cash->trap entry
gamma_day    = 0.2;              % Conditional daily removal probability
gamma_dt     = 1 - (1 - gamma_day)^(dt / 24); % Probability per time step

% Visualization
snapshot_days = [0, 2, 30];
snap_steps    = round(snapshot_days * 24 / dt) + 1;
plot_xlim     = [-20 120];       % Display window; does NOT confine movement
bin_width     = 2.5;             % Histogram width (m)
save_figure   = true;           % Set false to skip cropped PNG export

% Basic input checks
assert(k_cash > 0 && k_cash < 1, 'Require 0 < k_cash < 1.');
assert(rho > 0 && rho <= 1, 'Require 0 < rho <= 1.');
assert(abs(alpha) <= 1 && omega >= 0 && omega <= 1, ...
    'Invalid direction/persistence parameters.');
assert(beta_max >= 0 && beta_max <= 1 && nu >= 0, ...
    'Invalid spatial attraction parameters.');
assert(gamma_entry >= 0 && gamma_entry <= 1 && ...
       gamma_day >= 0 && gamma_day <= 1, 'Invalid removal probabilities.');
assert(all(snap_steps >= 1 & snap_steps <= p_steps), ...
    'All snapshots must fall within the simulation interval.');

%% 2. INITIAL CONDITIONS
% Uniformly distribute insects across the cash crop, x in [0, k*L].
% Give each insect an independent initial direction +/-1 with P=0.5.
x   = zone_left * rand(Npart, 1);
dir = sign(randn(Npart, 1));  % Independent, symmetric initial directions
dir(dir == 0) = 1;          % Extremely unlikely zero draw: choose +1

% Individual-level removal indicator. Once true, an insect never moves.
is_killed = false(Npart, 1);

time_hours = (0:p_steps-1) * dt;
time_days  = time_hours / 24;

%% 3. STORAGE AND INITIAL POPULATION ACCOUNTING
Nc_cash_pct     = zeros(1, p_steps);
Nc_trap_pct     = zeros(1, p_steps);
Nt_killed       = zeros(1, p_steps);
Nt_emigrated    = zeros(1, p_steps);
x_snap         = cell(1, numel(snap_steps));

Nc_cash_pct(1) = 100;            % All insects start in the cash crop
x_snap{1} = x;

%% 4. INDIVIDUAL-BASED SIMULATION
for a = 2:p_steps

    % Classify insects by their position at the START of the step.
    active      = ~is_killed;
    in_cash     = active & (x >= 0) & (x <= zone_left);
    in_trap     = active & (x > zone_left) & (x <= L);
    outside     = active & ((x < 0) | (x > L));
    was_in_cash = in_cash;       % Used to identify cash->trap entries

    % 4.1 Draw positive, independent half-normal movement step lengths.
    % Trap insects have a reduced movement scale, modelling arrestment.
    % Emigrated insects retain the same scale as insects in the cash crop.
    sigma_eff = sigma_cash * ones(Npart, 1);
    sigma_eff(in_trap) = sigma_trap;
    step_mag = abs(randn(Npart, 1)) .* sigma_eff;

    % 4.2 Construct spatial attraction for insects currently in cash crop.
    % No attraction is applied to insects in the trap or emigration zones.
    beta_eff = zeros(Npart, 1);
    beta_eff(in_cash) = beta_max * (x(in_cash) ./ zone_left).^nu;

    % 4.3 Update directions within the CASH CROP.
    % First generate persistence-driven directions:
    %   P(same heading) = (1+alpha)/2,
    %   P(reverse)      = (1-alpha)/2.
    dir_new = dir;
    keep_same = rand(Npart, 1) < (1 + alpha)/2;
    flip_idx = in_cash & ~keep_same;
    dir_new(flip_idx) = -dir_new(flip_idx);

    % Independently select which cash-crop insects follow the attraction
    % cue. For those insects, replace the persistence-driven direction by
    % a bias-driven draw: P(+1) = (1+beta(x))/2.
    do_bias = in_cash & (rand(Npart, 1) < omega);
    choose_toward = rand(Npart, 1) < (1 + beta_eff)/2;
    dir_new(do_bias) = 1;
    dir_new(do_bias & ~choose_toward) = -1;

    % 4.4 Update directions in the TRAP and EMIGRATION REGIONS.
    % Each insect in either region receives an independent +/-1 draw,
    % with probability 1/2 for each direction. The index vector maps
    % random draws to the corresponding individuals.
    trap_or_outside = in_trap | outside;
    idx = find(trap_or_outside);
    if ~isempty(idx)
        r = rand(numel(idx), 1);
        dir_new(idx) = 1;
        dir_new(idx(r > 0.5)) = -1;
    end
    dir_new(~active) = 0;

    % 4.5 Move every active insect by its signed step length.
    x(active) = x(active) + dir_new(active) .* step_mag(active);
    dir = dir_new;

    % 4.6 ENTRY-TRIGGERED REMOVAL: crossing the cash/trap interface
    % from the cash-crop side triggers a Bernoulli trial. The crossing
    % criterion also covers a rare step that overshoots the entire strip;
    % an insect cannot bypass the entry boundary by taking a long step.
    % A later cash-to-trap crossing generates another independent trial.
    entered_now = was_in_cash & (x > zone_left);
    kill_on_entry = entered_now & (rand(Npart, 1) < gamma_entry);
    is_killed(kill_on_entry) = true;

    % 4.7 WITHIN-TRAP REMOVAL: surviving insects in the trap are exposed
    % to the per-step equivalent of the daily removal probability.
    still_active = ~is_killed;
    in_trap_after_entry = still_active & (x > zone_left) & (x <= L);
    kill_inside = in_trap_after_entry & (rand(Npart, 1) < gamma_dt);
    is_killed(kill_inside) = true;

    % Exclude removed individuals permanently from positions/trajectories.
    just_killed = kill_on_entry | kill_inside;
    x(just_killed)   = NaN;
    dir(just_killed) = 0;

    % 4.8 Count mutually exclusive population states after both removals.
    alive = ~is_killed;
    Nc_cash_pct(a) = 100 * sum(alive & (x >= 0) & (x <= zone_left)) / Npart;
    Nc_trap_pct(a) = 100 * sum(alive & (x > zone_left) & (x <= L)) / Npart;
    Nt_killed(a)   = sum(is_killed);
    Nt_emigrated(a)= sum(alive & ((x < 0) | (x > L)));

    % Store positions for the specified day-2 and day-30 snapshots.
    snap_idx = find(snap_steps == a);
    if ~isempty(snap_idx)
        x_snap{snap_idx} = x;
    end
end

%% 5. CHECK POPULATION CONSERVATION AND REPORT RESULTS
Nt_killed_pct    = 100 * Nt_killed / Npart;
Nt_emigrated_pct = 100 * Nt_emigrated / Npart;

population_total = Nc_cash_pct + Nc_trap_pct + ...
                   Nt_killed_pct + Nt_emigrated_pct;
max_imbalance = max(abs(population_total - 100));
assert(max_imbalance < 1e-10, ...
    'Population accounting failed: four states do not sum to 100%.');

fprintf('Figure 2 simulation: N = %d, dt = %.2f hours, duration = %.0f days\n', ...
    Npart, dt, t_days);
fprintf('Max population-accounting discrepancy = %.3g percentage points\n', ...
    max_imbalance);
fprintf('Day 30: cash %.2f%% | trap %.2f%% | removed %.2f%% | emigrated %.2f%%\n', ...
    Nc_cash_pct(end), Nc_trap_pct(end), ...
    Nt_killed_pct(end), Nt_emigrated_pct(end));

%% 6. FIGURE: THREE SPATIAL HISTOGRAMS + POPULATION TIME SERIES
% Histograms display positions within plot_xlim only. Emigrants beyond
% this window remain in the simulated population and are counted in (d).

cash_col  = [0.20 0.52 0.20];    % Green
trap_col  = [0.53 0.25 0.67];    % Purple
emig_col  = [0.26 0.49 0.72];    % Blue
trap_shade= [0.86 0.78 0.88];   % Light purple
bin_edges = plot_xlim(1):bin_width:plot_xlim(2);

figure('Color', 'w', 'Position', [60 60 900 800]);
ymax = zeros(1, 3);

for panel = 1:3
    subplot(2, 2, panel);
    positions = x_snap{panel};
    positions = positions(~isnan(positions)); % Exclude permanently removed

    cash_data = positions(positions >= 0 & positions <= zone_left);
    trap_data = positions(positions > zone_left & positions <= L);
    emig_data = positions(positions < 0 | positions > L);

    histogram(cash_data, bin_edges, 'FaceColor', cash_col, ...
        'EdgeColor', 'none');
    hold on;
    histogram(trap_data, bin_edges, 'FaceColor', trap_col, ...
        'EdgeColor', 'none');
    histogram(emig_data, bin_edges, 'FaceColor', emig_col, ...
        'EdgeColor', 'none');

    y_current = ylim;
    ymax(panel) = y_current(2);
    title(sprintf('(%s) t = %g days', char('a' + panel-1), ...
        snapshot_days(panel)), 'FontName', 'Arial', 'FontWeight', 'normal');
    xlabel('Position x (m)', 'FontName', 'Arial');
    ylabel('Number of individuals', 'FontName', 'Arial');
    xlim(plot_xlim);
    set(gca, 'FontName', 'Arial', 'LineWidth', 1.2, 'Box', 'off');

    outside_plot = sum(positions < plot_xlim(1) | positions > plot_xlim(2));
    fprintf('Day %g: %d active insects beyond histogram window [%g,%g] m\n', ...
        snapshot_days(panel), outside_plot, plot_xlim(1), plot_xlim(2));
    hold off;
end

% Use a common y-axis across the three histogram panels and indicate
% the trap-zone extent and the field boundaries.
common_ylim = [0 max(ymax)];
for panel = 1:3
    subplot(2, 2, panel);
    ylim(common_ylim);
    hold on;

    patch([zone_left L L zone_left], ...
          [common_ylim(1) common_ylim(1) common_ylim(2) common_ylim(2)], ...
          trap_shade, 'EdgeColor', 'none', 'FaceAlpha', 0.25);

    plot([0 0], common_ylim, 'k--', 'LineWidth', 1);
    plot([zone_left zone_left], common_ylim, 'k--', 'LineWidth', 1);
    plot([L L], common_ylim, 'k--', 'LineWidth', 1);
    hold off;
end

% Panel (d): percentages of the initial Npart in each state through time.
subplot(2, 2, 4);
axis square;
hold on;
plot(time_days, Nc_cash_pct, 'Color', cash_col, 'LineWidth', 1.5);
plot(time_days, Nc_trap_pct, 'Color', trap_col, 'LineWidth', 1.5);
plot(time_days, Nt_killed_pct, 'Color', [0.85 0 0], 'LineWidth', 1.5);
plot(time_days, Nt_emigrated_pct, 'Color', emig_col, 'LineWidth', 1.5);

xlabel('Time (days)', 'FontName', 'Arial');
ylabel('Population (%)', 'FontName', 'Arial');
legend({'Cash crop', 'Trap crop', 'Removed', 'Emigrated'}, ...
    'Location', 'northeast', 'Box', 'off');
xlim([0 t_days]);
ylim([0 100]);
set(gca, 'FontName', 'Arial', 'LineWidth', 1.5, 'Box', 'off');
hold off;

% Tight PNG export: trim only the outside white margin; retain all
% axis labels, legends, tick labels and colour bars. No plotted data change.
if save_figure
    png_file = 'Fig_2.png';
    set(gcf, 'PaperPositionMode', 'auto');  % Preserve on-screen layout
    print(gcf, png_file, '-dpng', '-r300');

    % Crop the exported pixels, not the MATLAB figure or the axes.
    % This uses base MATLAB imread/imwrite; no toolbox is required.
    rgb = imread(png_file);
    nonwhite = any(rgb(:,:,1:min(3,size(rgb,3))) < 250, 3);
    [rows, cols] = find(nonwhite);
    if ~isempty(rows)
        padding = 18;  % Small, consistent safety margin (pixels)
        top    = max(1, min(rows) - padding);
        bottom = min(size(rgb,1), max(rows) + padding);
        left   = max(1, min(cols) - padding);
        right  = min(size(rgb,2), max(cols) + padding);
        imwrite(rgb(top:bottom,left:right,:), png_file, 'png');
    end
    fprintf('Saved tightly cropped figure: %s\n', png_file);
end
