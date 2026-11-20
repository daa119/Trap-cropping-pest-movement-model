%% FIGURE 6 | INTERCEPTION, REMOVAL AND ARRESTMENT DYNAMICS
% Insect movement and spatial dynamics determine trap crop effectiveness
% Ecological Modelling (2026)
%
% PURPOSE
%   Compare six scenarios that vary entry-triggered interception,
%   mortality during trap-crop residence and behavioural arrestment.
%   Simulate Npart insects for 30 days and plot four population states:
%       (a) Cash-crop occupancy
%       (b) Trap-crop occupancy
%       (c) Cumulative removal
%       (d) Current emigration
%
% MODEL GEOMETRY AND POPULATION STATES
%   A one-dimensional field of length L is divided into a cash crop
%   (0 <= x <= k*L) and a terminal trap crop (k*L < x <= L).
%   Insects outside the field (x < 0 or x > L) remain active and can
%   return. Removed insects stop moving permanently. At every time point,
%   the four mutually exclusive states account for all Npart insects.
%
% MOVEMENT RULES
%   Position is updated by x_i(t) = x_i(t-1) + delta_i(t)*s_i(t).
%   Positive step lengths are half-normal:
%       s_i = sigma_eff * abs(Z),  Z ~ Normal(0,1).
%   In the cash crop, movement direction is chosen by either:
%     (1) persistence, with probability (1-omega), for which
%         P(same direction) = (1+alpha)/2, or
%     (2) attraction, with probability omega, for which
%         beta(x) = beta_max*(x/(k*L))^nu and
%         P(rightward) = (1+beta(x))/2.
%   Inside the trap crop and outside the field, reorientation is unbiased:
%       P(rightward) = P(leftward) = 1/2 independently at each step.
%   The step-length scale in the trap crop is rho*sigma_cash; elsewhere
%   it is sigma_cash. Thus rho=1 means no arrestment.
%
% REMOVAL RULES
%   An insect crossing from the cash crop into the trap crop is removed
%   with probability gamma_entry. A rare step spanning the full trap crop
%   still crosses this interface and is eligible for entry removal.
%   Survivors currently within the trap crop are also subject to a
%   per-step removal probability
%       gamma_dt = 1 - (1-gamma_day)^(dt/24),  dt measured in hours.
%   Removal on entry and during residence are separate mechanisms.
%
% SIMULATION DESIGN
%   Each scenario uses one stochastic realization of 5000 insects, rather
%   than an average over independent simulations. A fixed random seed
%   ensures reproducibility on a compatible MATLAB installation. The
%   baseline scenario uses the same parameters and random-draw sequence
%   as the corrected Figure 2 script.
%
% COMPATIBILITY
%   MATLAB R2016a or later; no additional toolboxes are required.

clear; clc; close all;
rng(1, 'twister');

%% 1. FIXED PARAMETERS (BASELINE AS IN FIGURE 2)
dt          = 0.5;                  % Time step (hours)
Npart       = 5000;                 % Initial number of insects
t_days      = 30;                   % Total duration (days)
p_steps     = round(t_days*24/dt)+1; % Include initial time point
sigma_cash  = 1;                    % Half-normal step-length scale (m)
L           = 100;                  % Total field length (m)
k_cash      = 0.8;                  % Fraction of field occupied by cash crop
zone_left   = k_cash*L;            % Cash/trap boundary (m)

omega       = 0.5;                  % Probability of using environmental bias
beta_max    = 0.1;                  % Maximum attraction strength
nu          = 1;                    % Linear spatial attraction profile
alpha       = 0.2;                  % Directional persistence

time_hours  = (0:p_steps-1)*dt;
time_days   = time_hours/24;

%% 2. SCENARIOS (ENTRY REMOVAL, WITHIN-TRAP REMOVAL, ARRESTMENT)
% Columns correspond to the six scenarios listed in the Figure 6 caption.
gamma_entry_list = [0.1, 0,   0.1, 0.1, 0.1, 1];
gamma_day_list   = [0.2, 0.2, 0,   0.2, 0.01, 0.2];
rho_list         = [0.1, 0.1, 0.1, 1,   0.01, 0.1];

labels = { ...
    'Baseline', ...
    'No entry interception', ...
    'No within-trap crop removal', ...
    'No arrestment', ...
    'Strong arrestment and weak lethality', ...
    'Perfect entry interception'};

nCases = numel(labels);

% Verify the parameter ranges used in every scenario.
assert(dt > 0 && Npart > 0 && L > 0 && sigma_cash > 0, ...
    'Time, population, field length and movement scale must be positive.');
assert(k_cash > 0 && k_cash < 1, 'Require 0 < k_cash < 1.');
assert(abs(alpha) <= 1 && omega >= 0 && omega <= 1, ...
    'Invalid persistence or environmental bias weight.');
assert(beta_max >= 0 && beta_max <= 1 && nu >= 0, ...
    'Invalid spatial attraction parameters.');
assert(all(gamma_entry_list >= 0 & gamma_entry_list <= 1) && ...
       all(gamma_day_list >= 0 & gamma_day_list <= 1), ...
    'Removal probabilities must be between zero and one.');
assert(all(rho_list > 0 & rho_list <= 1), 'Require 0 < rho <= 1.');
assert(nCases == numel(gamma_entry_list) && ...
       nCases == numel(gamma_day_list) && nCases == numel(rho_list), ...
       'Every scenario needs a complete parameter specification.');

%% 3. TIME-SERIES STORAGE
cash_all = zeros(nCases, p_steps);
trap_all = zeros(nCases, p_steps);
kill_all = zeros(nCases, p_steps);
emig_all = zeros(nCases, p_steps);

%% 4. SIMULATE EACH SCENARIO
for icase = 1:nCases
    gamma_entry = gamma_entry_list(icase);
    gamma_day   = gamma_day_list(icase);
    rho         = rho_list(icase);

    sigma_trap  = rho*sigma_cash;
    gamma_dt    = 1 - (1-gamma_day)^(dt/24);

    % Insects initially occupy the cash crop uniformly, with independent
    % and equally likely leftward and rightward headings.
    x = zone_left*rand(Npart,1);
    dir = sign(randn(Npart,1));
    dir(dir == 0) = 1;
    is_killed = false(Npart,1);

    % All insects are initially alive and inside the cash crop.
    cash_all(icase,1) = 100;
    trap_all(icase,1) = 0;
    kill_all(icase,1) = 0;
    emig_all(icase,1) = 0;

    for a = 2:p_steps
        active = ~is_killed;

        % Classify by position at the START of this movement step.
        in_cash = active & (x >= 0) & (x <= zone_left);
        in_trap = active & (x > zone_left) & (x <= L);
        outside = active & ((x < 0) | (x > L));
        was_in_cash = in_cash;

        % Draw positive half-normal steps. Arrestment applies only inside
        % the trap crop, not in the cash crop or emigration regions.
        sigma_eff = sigma_cash*ones(Npart,1);
        sigma_eff(in_trap) = sigma_trap;
        step_mag = abs(randn(Npart,1)).*sigma_eff;

        % Attraction varies with position within the cash crop and is
        % switched off everywhere else.
        beta_eff = zeros(Npart,1);
        beta_eff(in_cash) = beta_max*(x(in_cash)./zone_left).^nu;

        % Persistence-driven headings for insects in the cash crop.
        dir_new = dir;
        keep_same = rand(Npart,1) < (1+alpha)/2;
        flip_idx = in_cash & ~keep_same;
        dir_new(flip_idx) = -dir_new(flip_idx);

        % Independently choose which cash-crop insects follow the cue.
        do_bias = in_cash & (rand(Npart,1) < omega);
        choose_toward = rand(Npart,1) < (1+beta_eff)/2;
        dir_new(do_bias) = 1;
        dir_new(do_bias & ~choose_toward) = -1;

        % Unbiased simple random walk in trap and emigration regions.
        % Use explicit indices so each coin toss updates the correct insect.
        idx = find(in_trap | outside);
        if ~isempty(idx)
            r = rand(numel(idx),1);
            dir_new(idx) = 1;
            dir_new(idx(r > 0.5)) = -1;
        end
        dir_new(~active) = 0;

        % Move all active insects once.
        x(active) = x(active) + dir_new(active).*step_mag(active);
        dir = dir_new;

        % Entry-triggered removal whenever a step crosses the cash/trap
        % boundary from the cash crop, even if the step ends past x=L.
        entered_now = was_in_cash & (x > zone_left);
        kill_on_entry = entered_now & (rand(Npart,1) < gamma_entry);
        is_killed(kill_on_entry) = true;

        % Within-trap removal among surviving insects CURRENTLY in trap.
        still_active = ~is_killed;
        in_trap_after_entry = still_active & (x > zone_left) & (x <= L);
        kill_inside = in_trap_after_entry & ...
                      (rand(Npart,1) < gamma_dt);
        is_killed(kill_inside) = true;

        % Permanently remove killed insects from subsequent movements.
        just_killed = kill_on_entry | kill_inside;
        x(just_killed) = NaN;
        dir(just_killed) = 0;

        % Count four non-overlapping states at the END of this step.
        alive = ~is_killed;
        nCash = sum(alive & (x >= 0) & (x <= zone_left));
        nTrap = sum(alive & (x > zone_left) & (x <= L));
        nRemoved = sum(is_killed);
        nEmigrated = sum(alive & ((x < 0) | (x > L)));

        % Every original insect must be in exactly one state.
        assert(nCash+nTrap+nRemoved+nEmigrated == Npart, ...
            'Population accounting failed at scenario %d, step %d.', ...
            icase, a);

        cash_all(icase,a) = 100*nCash/Npart;
        trap_all(icase,a) = 100*nTrap/Npart;
        kill_all(icase,a) = 100*nRemoved/Npart;
        emig_all(icase,a) = 100*nEmigrated/Npart;
    end

    fprintf('%-38s | day 30: cash %5.1f%%, trap %5.1f%%, removed %5.1f%%, emigrated %5.1f%%\n', ...
        labels{icase}, cash_all(icase,end), trap_all(icase,end), ...
        kill_all(icase,end), emig_all(icase,end));
end

%% 5. VERIFY TOTAL POPULATION AT EVERY RECORDED TIME
population_totals = cash_all + trap_all + kill_all + emig_all;
max_imbalance = max(abs(population_totals(:)-100));
assert(max_imbalance < 1e-10, ...
    'The four population states must add to 100%% at every time point.');
fprintf('Maximum population-accounting discrepancy: %.3g percentage points\n', ...
    max_imbalance);

%% 6. PLOT THE FOUR POPULATION STATES
% Identical line colours and styles identify the six scenarios in every
% panel. The legend is placed in panel (a), as in the manuscript figure.
colors = [
    0.00 0.00 0.00
    0.90 0.62 0.00
    0.00 0.45 0.70
    0.00 0.62 0.45
    0.80 0.47 0.65
    0.84 0.37 0.00
];
lineStyles = {'-','--',':','-.','--','-'};

% Consistent 2-by-2 plot-box dimensions and spacing, matched to Figure 2.
% Units are pixels; all four axes have the same width and height.
% This is strictly a plotting adjustment; the simulation is unchanged.
panel_positions_px = [ ...
     70 460 335 290; ...  % (a) upper left
    535 460 335 290; ...  % (b) upper right
     70  65 335 290; ...  % (c) lower left
    535  65 335 290];     % (d) lower right

figure('Color','w','Position',[60 60 900 800]);

axCash = subplot(2,2,1);
hLegend = gobjects(nCases,1);
hold(axCash,'on');
for i = 1:nCases
    hLegend(i) = plot(time_days, cash_all(i,:), ...
        'Color',colors(i,:), 'LineStyle',lineStyles{i},'LineWidth',2);
end
ylabel('Cash-crop occupancy (%)');
title('(a) Cash-crop occupancy','FontWeight','normal');
hold(axCash,'off');

axTrap = subplot(2,2,2);
hold(axTrap,'on');
for i = 1:nCases
    plot(time_days, trap_all(i,:), ...
        'Color',colors(i,:), 'LineStyle',lineStyles{i},'LineWidth',2);
end
ylabel('Trap-crop occupancy (%)');
title('(b) Trap-crop occupancy','FontWeight','normal');
hold(axTrap,'off');

axRemoved = subplot(2,2,3);
hold(axRemoved,'on');
for i = 1:nCases
    plot(time_days, kill_all(i,:), ...
        'Color',colors(i,:), 'LineStyle',lineStyles{i},'LineWidth',2);
end
xlabel('Time (days)');
ylabel('Removed (%)');
title('(c) Removed','FontWeight','normal');
hold(axRemoved,'off');

axEmigrated = subplot(2,2,4);
hold(axEmigrated,'on');
for i = 1:nCases
    plot(time_days, emig_all(i,:), ...
        'Color',colors(i,:), 'LineStyle',lineStyles{i},'LineWidth',2);
end
xlabel('Time (days)');
ylabel('Emigrated (%)');
title('(d) Emigrated','FontWeight','normal');
hold(axEmigrated,'off');

legend(axCash,hLegend,labels, ...
    'Location','northeast','Box','off', ...
    'FontName','Arial','FontSize',8);


% Lock the plotting boxes to the common four-panel layout after legends
% (and colour bars, if present) have been created.
panel_axes = [axCash axTrap axRemoved axEmigrated];
for iPanel = 1:4
    set(panel_axes(iPanel), 'Units', 'pixels', ...
        'Position', panel_positions_px(iPanel,:));
end

fprintf('\nFigure 6 simulation complete.\n');
%% AUTOMATIC TIGHT PNG EXPORT (VISUAL OUTPUT ONLY)
% Removes blank canvas outside the plotted content, including labels,
% legend and colour bar. Does not change simulations or axes.
save_figure = true;    % Set false to display without exporting
if save_figure
    png_file = 'Fig_6.png';
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
