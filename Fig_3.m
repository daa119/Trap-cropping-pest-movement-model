%% FIGURE 3 | SPATIAL ATTRACTION AND TRAP-CROP ALLOCATION
% Insect movement and spatial dynamics determine trap crop effectiveness
% Ecological Modelling (2026)
%
% PURPOSE
%   Evaluate the fractions of insects in four population states after
%   15 days as functions of maximum attraction beta_max, attraction-profile
%   exponent nu, and the fraction k of the field planted as cash crop.
%   The script generates the four panels of Figure 3:
%       (a) Cash-crop occupancy
%       (b) Trap-crop occupancy
%       (c) Cumulative removal
%       (d) Current emigration
%
% MODEL GEOMETRY AND STATES
%   Field length L, with a cash crop in 0 <= x <= k*L and a trap crop in
%   k*L < x <= L. Insects outside [0,L] are temporarily emigrated and may
%   return. Removed insects never move again. These four population states
%   are mutually exclusive and sum to the initial population Npart.
%
% MOVEMENT MODEL
%   Individual movement uses x_i(t+1) = x_i(t) + delta_i(t+1)*s_i(t+1).
%   Positive step lengths are drawn from a half-normal distribution:
%       s = sigma_eff * abs(Z),   Z ~ Normal(0,1).
%   In the cash crop, a new direction is selected from either:
%       persistence, with probability 1-omega, or
%       attraction bias, with probability omega.
%   For persistence, P(same direction) = (1+alpha)/2.
%   For attraction,  P(rightward) = (1+beta(x))/2, where
%       beta(x) = beta_max * (x/(k*L))^nu,  0 <= x <= k*L.
%   In the trap crop and beyond the field, each direction is chosen
%   independently, with P(leftward) = P(rightward) = 1/2.
%   Trap-crop movement uses the reduced step scale rho*sigma_cash.
%
% REMOVAL MODEL
%   On crossing from the cash crop into the trap crop, an insect is removed
%   with probability gamma_entry. Surviving insects inside the trap crop
%   are then subject to per-step removal with probability
%       gamma_dt = 1 - (1-gamma_day)^(dt/24), where dt is in hours.
%   Entry removal can occur again if an insect returns to the cash crop
%   and subsequently enters the trap crop a second time.
%
% SIMULATION DESIGN
%   One simulation of Npart independent insects is run for every combination
%   of k, nu and beta_max. This is a stochastic parameter sweep, not an
%   average over replicate simulations. A fixed random seed makes the
%   script reproducible on a compatible MATLAB installation.
%
% COMPATIBILITY
%   MATLAB R2016a or later; no additional toolboxes required.

clear; clc; close all;
rng(1, 'twister');

%% 1. FIXED MODEL PARAMETERS (AS IN FIGURE 2)
dt          = 0.5;                  % Duration of a movement step (hours)
Npart       = 5000;                 % Initial number of insects
t_days      = 15;                   % Evaluation time (days)
p_steps     = round(t_days*24/dt)+1; % Include the initial time point
L           = 100;                  % Total field length (m)
sigma_cash  = 1;                    % Half-normal step scale in cash crop (m)

alpha       = 0.2;                  % Persistence correlation parameter
omega       = 0.5;                  % Probability of a cue-driven direction
rho         = 0.1;                  % Step-length reduction in trap crop
sigma_trap  = rho*sigma_cash;       % Trap-crop step-length scale (m)

gamma_entry = 0.1;                  % Removal probability on cash-to-trap entry
gamma_day   = 0.2;                  % Removal probability for 24 h in trap
gamma_dt    = 1-(1-gamma_day)^(dt/24); % Per-step removal probability

%% 2. PARAMETERS VARIED IN FIGURE 3
k_vals      = [0.7, 0.9];           % Cash-crop fractions: 30% and 10% trap
nu_vals     = [0, 1, 2];            % Uniform, linear and concentrated bias
beta_vals   = 0:0.05:1;             % Maximum attraction strengths

nK = numel(k_vals);
nV = numel(nu_vals);
nB = numel(beta_vals);

% Store percentages as arrays indexed by (k, nu, beta_max).
final_cash_pct      = zeros(nK, nV, nB);
final_trap_pct      = zeros(nK, nV, nB);
final_removed_pct   = zeros(nK, nV, nB);
final_emigrated_pct = zeros(nK, nV, nB);

% Fundamental checks of the model's admissible parameter ranges.
assert(dt > 0 && Npart > 0 && L > 0 && sigma_cash > 0, ...
    'Time step, insect number, length and movement scale must be positive.');
assert(all(k_vals > 0 & k_vals < 1), 'Require 0 < k < 1.');
assert(all(nu_vals >= 0), 'Require nu >= 0.');
assert(all(beta_vals >= 0 & beta_vals <= 1), ...
    'Require 0 <= beta_max <= 1.');
assert(abs(alpha) <= 1 && omega >= 0 && omega <= 1, ...
    'Invalid alpha or omega.');
assert(rho > 0 && rho <= 1, 'Require 0 < rho <= 1.');
assert(gamma_entry >= 0 && gamma_entry <= 1 && ...
       gamma_day >= 0 && gamma_day <= 1, ...
       'Removal probabilities must lie between 0 and 1.');

%% 3. STOCHASTIC PARAMETER SWEEP
% The seed is fixed once before the sweep. Each parameter combination
% therefore receives its own sequence of random movement and removal draws.
max_population_imbalance = 0;

for iK = 1:nK
    k_cash = k_vals(iK);
    zone_left = k_cash*L;           % Cash-crop / trap-crop boundary
    fprintf('Simulating cash-crop fraction k = %.2f ...\n', k_cash);

    for iV = 1:nV
        nu = nu_vals(iV);

        for iB = 1:nB
            beta_max = beta_vals(iB);

            % 3.1 Initialise the insects in the cash crop.
            x   = zone_left*rand(Npart,1);
            dir = sign(randn(Npart,1));
            dir(dir == 0) = 1;      % Assign +1 for the extremely rare zero
            is_killed = false(Npart,1);

            % 3.2 Perform the required number of movement steps.
            for a = 2:p_steps
                active = ~is_killed;

                % Spatial states at the beginning of the step.
                in_cash = active & (x >= 0) & (x <= zone_left);
                in_trap = active & (x > zone_left) & (x <= L);
                outside = active & ((x < 0) | (x > L));
                was_in_cash = in_cash;

                % Half-normal movement steps; reduced scale only in trap.
                sigma_eff = sigma_cash*ones(Npart,1);
                sigma_eff(in_trap) = sigma_trap;
                step_mag = abs(randn(Npart,1)).*sigma_eff;

                % Spatially varying attraction within the cash crop.
                beta_eff = zeros(Npart,1);
                beta_eff(in_cash) = beta_max * ...
                    (x(in_cash)./zone_left).^nu;

                % Persistence-based directions in the cash crop.
                dir_new = dir;
                keep_same = rand(Npart,1) < (1+alpha)/2;
                reverse = in_cash & ~keep_same;
                dir_new(reverse) = -dir_new(reverse);

                % Independently select which cash-crop insects use cue bias.
                % Bias overrides persistence for the selected insects.
                use_bias = in_cash & (rand(Npart,1) < omega);
                choose_right = rand(Npart,1) < (1+beta_eff)/2;
                dir_new(use_bias) = 1;
                dir_new(use_bias & ~choose_right) = -1;

                % Unbiased reorientation in trap and emigration regions.
                idx = find(in_trap | outside);
                if ~isempty(idx)
                    r = rand(numel(idx),1);
                    dir_new(idx) = 1;
                    dir_new(idx(r > 0.5)) = -1;
                end
                dir_new(~active) = 0;

                % Update the positions of all active insects.
                x(active) = x(active) + ...
                    dir_new(active).*step_mag(active);
                dir = dir_new;

                % Cash-to-trap boundary crossings trigger entry removal.
                % The condition also covers a step crossing the full strip
                % in a single move (rare at the specified movement scales).
                entered_now = was_in_cash & (x > zone_left);
                kill_on_entry = entered_now & ...
                    (rand(Npart,1) < gamma_entry);
                is_killed(kill_on_entry) = true;

                % Survivors currently inside the trap face per-step removal.
                still_active = ~is_killed;
                in_trap_after_entry = still_active & ...
                    (x > zone_left) & (x <= L);
                kill_inside = in_trap_after_entry & ...
                    (rand(Npart,1) < gamma_dt);
                is_killed(kill_inside) = true;

                % Prevent subsequent movement of permanently removed insects.
                just_killed = kill_on_entry | kill_inside;
                x(just_killed) = NaN;
                dir(just_killed) = 0;
            end

            % 3.3 Classify all insects at the final evaluation time.
            alive = ~is_killed;
            nCash = sum(alive & (x >= 0) & (x <= zone_left));
            nTrap = sum(alive & (x > zone_left) & (x <= L));
            nRemoved = sum(is_killed);
            nEmigrated = sum(alive & ((x < 0) | (x > L)));

            % Four exclusive states must account for every initial insect.
            imbalance = abs(nCash+nTrap+nRemoved+nEmigrated-Npart);
            max_population_imbalance = max(max_population_imbalance, imbalance);
            assert(imbalance == 0, 'Population states fail to partition Npart.');

            final_cash_pct(iK,iV,iB) = 100*nCash/Npart;
            final_trap_pct(iK,iV,iB) = 100*nTrap/Npart;
            final_removed_pct(iK,iV,iB) = 100*nRemoved/Npart;
            final_emigrated_pct(iK,iV,iB) = 100*nEmigrated/Npart;
        end
    end
end

fprintf('Simulation complete: %d scenarios.\n', nK*nV*nB);
fprintf('Maximum population-count imbalance: %g insects.\n', ...
    max_population_imbalance);

%% 4. FOUR-PANEL FIGURE
% Solid, dash-dot and dotted lines denote nu = 0, 1 and 2, respectively.
% Light/dark shades denote k = 0.7 / k = 0.9 in each panel.
% Panel colours denote population state (cash/trap/removed/emigrated).

metric_colors = {
    {[0.55 0.78 0.55], [0.20 0.52 0.20]}, ... % Green: cash
    {[0.73 0.58 0.80], [0.53 0.25 0.67]}, ... % Purple: trap
    {[0.95 0.55 0.55], [0.85 0.00 0.00]}, ... % Red: removed
    {[0.56 0.72 0.88], [0.26 0.49 0.72]}  ... % Blue: emigrated
};
line_styles = {'-', '-.', ':'};
metric_names = {'Cash crop occupancy', 'Trap crop occupancy', ...
                'Removed', 'Emigrated'};
data_arrays = {final_cash_pct, final_trap_pct, ...
               final_removed_pct, final_emigrated_pct};
reference_ylim = [100 20 100 20]; % Fixed panel y-axis maxima (a)-(d)
save_figure = true;              % Set false to skip cropped PNG export

% Consistent 2-by-2 plot-box dimensions and spacing, matched to Figure 2.
% Units are pixels; all four axes have the same width and height.
% This is strictly a plotting adjustment; the simulation is unchanged.
panel_positions_px = [ ...
     70 460 335 290; ...  % (a) upper left
    535 460 335 290; ...  % (b) upper right
     70  65 335 290; ...  % (c) lower left
    535  65 335 290];     % (d) lower right

figure('Color','w','Position',[60 60 900 800]);
legend_handles = gobjects(nV,1);
panel_axes = gobjects(4,1);

for iMetric = 1:4
    panel_axes(iMetric) = subplot(2,2,iMetric);
    hold on;

    results = data_arrays{iMetric};
    shades = metric_colors{iMetric};

    for iK = 1:nK
        for iV = 1:nV
            y = squeeze(results(iK,iV,:));
            plot(beta_vals,y, ...
                'LineStyle',line_styles{iV}, ...
                'LineWidth',2, ...
                'Color',shades{iK}, ...
                'HandleVisibility','off');
        end
    end

    % Create black line samples for the three nu legend entries.
    if iMetric == 1
        for iV = 1:nV
            legend_handles(iV) = plot(nan,nan, ...
                'LineStyle',line_styles{iV}, ...
                'LineWidth',2, 'Color','k');
        end
    end

    if iMetric >= 3
        xlabel('\beta_{max} (maximum attraction strength)', ...
            'FontName','Arial','FontSize',12,'Interpreter','tex');
    end
    if iMetric == 1 || iMetric == 3
        ylabel('Proportion of initial population (%)', ...
            'FontName','Arial','FontSize',12);
    end

    title(sprintf('(%s) %s',char('a'+iMetric-1),metric_names{iMetric}), ...
        'FontName','Arial','FontSize',12,'FontWeight','normal');
    xlim([0 1]);

    % Fixed y-axis limits for this panel.
    ylim([0 reference_ylim(iMetric)]);

    grid off;
    set(gca,'FontName','Arial','FontSize',11,'Box','off', ...
        'TickDir','out','LineWidth',1.5, ...
        'XMinorTick','off','YMinorTick','off');
    hold off;
end

legend_labels = cell(1,nV);
for iV = 1:nV
    legend_labels{iV} = sprintf('\\nu = %d',nu_vals(iV));
end
legend(legend_handles,legend_labels,'Interpreter','tex', ...
       'FontName','Arial','FontSize',10,'Box','off', ...
       'Location','best');


% Lock the plotting boxes to the common four-panel layout after legends
% (and colour bars, if present) have been created.
for iPanel = 1:4
    set(panel_axes(iPanel), 'Units', 'pixels', ...
        'Position', panel_positions_px(iPanel,:));
end

% Tight PNG export: trim only the outside white margin; retain all
% axis labels, legends, tick labels and colour bars. No plotted data change.
if save_figure
    png_file = 'Fig_3.png';
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
