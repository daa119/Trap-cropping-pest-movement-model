%% FIGURE 7 | EFFECT OF TRAP-CROP WIDTH ON POPULATION OUTCOMES
% Insect movement and spatial dynamics determine trap crop effectiveness
% Ecological Modelling (2026)
%
% PURPOSE
%   Evaluate how the fraction of the field assigned to the cash crop, k,
%   affects the four population states at t = 15 days under six scenarios
%   with different entry interception, within-trap removal and arrestment.
%   The trap crop occupies the terminal fraction (1-k) of a field of
%   fixed total length L. Each (scenario,k) pair is simulated separately.
%
% GEOMETRY AND POPULATION STATES
%   Cash crop:         0 <= x <= k*L
%   Trap crop:         k*L < x <= L
%   Emigration zones:  x < 0 or x > L
%   Removed insects are inactive; emigrated insects remain active and
%   can subsequently re-enter the field. At any time, each original
%   insect belongs to exactly one of these four states.
%
% INDIVIDUAL MOVEMENT
%   The signed step update is x_i(t) = x_i(t-1) + delta_i(t)*s_i(t).
%   Positive step lengths are half-normal:
%       s_i = sigma_region*abs(Z),  Z ~ Normal(0,1).
%   Within the cash crop, select one of two direction mechanisms:
%     - Directional persistence with probability (1-omega):
%         P(same direction) = (1+alpha)/2.
%     - Attraction with probability omega:
%         beta(x) = beta_max*(x/(k*L))^nu,
%         P(rightward) = (1+beta(x))/2.
%   In the trap crop AND in both emigration zones, directions follow an
%   unbiased simple random walk: P(leftward) = P(rightward) = 1/2.
%   In the trap crop only, movement scale is reduced to rho*sigma_cash.
%
% REMOVAL MECHANISMS
%   At every cash-to-trap interface crossing, remove with probability
%   gamma_entry. This includes a rare step that spans the full trap crop
%   and ends beyond x=L. Entry interception does not apply to crossings
%   from the right-hand emigration region into the trap crop.
%   During trap-crop residence, surviving insects are independently
%   removed with probability gamma_dt per time step, where
%       gamma_dt = 1-(1-gamma_day)^(dt/24).
%
% DESIGN AND REPRODUCIBILITY
%   Six scenarios use the same definitions and baseline parameters as
%   Figure 6, but the endpoint here is 15 days and k varies from 0.70 to
%   0.95. There is one stochastic run of Npart=5000 insects per parameter
%   combination, not an average across independent runs. The random seed
%   makes this sweep reproducible. Because random draws are consumed in a
%   different order from Figure 6, matching parameters need not yield
%   identical numerical percentages across the two separate scripts.
%
% COMPATIBILITY
%   MATLAB R2016a or later, without additional toolboxes.

clear; clc; close all;
rng(1, 'twister');

%% 1. FIXED MODEL PARAMETERS
Npart      = 5000;                 % Number of insects initially in cash crop
L          = 100;                  % Total field length (m)
dt         = 0.5;                  % Movement time step (hours)
t_days     = 15;                   % Time of population-state snapshot (days)
p_steps    = round(t_days*24/dt)+1; % Include the initial condition
sigma_cash = 1;                    % Half-normal scale in cash/outside (m)

alpha      = 0.2;                  % Directional persistence in cash crop
omega      = 0.5;                  % Frequency of choosing attraction in cash
beta_max   = 0.1;                  % Attraction magnitude at x=k*L
nu         = 1;                    % Linear spatial attraction profile

%% 2. CASH-CROP FRACTION AND SCENARIO SWEEPS
% k is the CASH-crop fraction; trap-crop width is L*(1-k).
k_vals = 0.70:0.05:0.95;
nK     = numel(k_vals);

% The six parameter combinations match those in Figure 6.
gamma_entry_list = [0.1, 0,   0.1, 0.1, 0.1, 1];
gamma_day_list   = [0.2, 0.2, 0,   0.2, 0.01, 0.2];
rho_list         = [0.1, 0.1, 0.1, 1,   0.01, 0.1];

labels = { ...
    'Baseline', ...
    'No entry interception', ...
    'No within-trap removal', ...
    'No arrestment', ...
    'Strong arrestment and weak lethality', ...
    'Perfect entry interception'};

nCases = numel(labels);

% Check the admissible parameter ranges before simulation.
assert(Npart > 0 && L > 0 && dt > 0 && sigma_cash > 0, ...
    'Population, geometry, time step and motility must be positive.');
assert(all(k_vals > 0 & k_vals < 1), 'Require 0 < k < 1.');
assert(abs(alpha) <= 1 && omega >= 0 && omega <= 1, ...
    'Invalid persistence or environmental bias weight.');
assert(beta_max >= 0 && beta_max <= 1 && nu >= 0, ...
    'Invalid attraction parameters.');
assert(all(gamma_entry_list >= 0 & gamma_entry_list <= 1) && ...
       all(gamma_day_list >= 0 & gamma_day_list <= 1), ...
    'Removal probabilities must lie between 0 and 1.');
assert(all(rho_list > 0 & rho_list <= 1), 'Require 0 < rho <= 1.');
assert(nCases == numel(gamma_entry_list) && ...
       nCases == numel(gamma_day_list) && ...
       nCases == numel(rho_list), ...
    'Each scenario requires all three mechanism parameters.');

%% 3. FINAL POPULATION PERCENTAGES
cash_final = zeros(nCases,nK);
trap_final = zeros(nCases,nK);
kill_final = zeros(nCases,nK);
emig_final = zeros(nCases,nK);

%% 4. SIMULATE EACH SCENARIO AND TRAP-CROP WIDTH
for icase = 1:nCases
    gamma_entry = gamma_entry_list(icase);
    gamma_day   = gamma_day_list(icase);
    rho         = rho_list(icase);

    sigma_trap  = rho*sigma_cash;
    gamma_dt    = 1-(1-gamma_day)^(dt/24);

    for ik = 1:nK
        k_cash    = k_vals(ik);
        zone_left = k_cash*L;   % Interface between cash and trap crops

        % Uniform initial positions across the CASH crop for this k.
        x   = zone_left*rand(Npart,1);
        dir = sign(randn(Npart,1));
        dir(dir == 0) = 1;
        is_killed = false(Npart,1);

        for a = 2:p_steps
            active = ~is_killed;

            % Determine region at the START of the step. These masks are
            % disjoint and exclude permanently removed insects.
            in_cash = active & (x >= 0) & (x <= zone_left);
            in_trap = active & (x > zone_left) & (x <= L);
            outside = active & ((x < 0) | (x > L));
            was_in_cash = in_cash;

            % Positive half-normal steps; arrestment only within trap.
            sigma_eff = sigma_cash*ones(Npart,1);
            sigma_eff(in_trap) = sigma_trap;
            step_mag = abs(randn(Npart,1)).*sigma_eff;

            % Spatial attraction is defined only for insects in CASH.
            beta_eff = zeros(Npart,1);
            beta_eff(in_cash) = ...
                beta_max*(x(in_cash)./zone_left).^nu;

            % Persistence-based direction for cash-crop insects.
            dir_new = dir;
            keep_same = rand(Npart,1) < (1+alpha)/2;
            flip_idx = in_cash & ~keep_same;
            dir_new(flip_idx) = -dir_new(flip_idx);

            % Independently select which cash-crop individuals follow
            % the attraction cue rather than the persistence component.
            do_bias = in_cash & (rand(Npart,1) < omega);
            choose_toward = rand(Npart,1) < (1+beta_eff)/2;
            dir_new(do_bias) = 1;
            dir_new(do_bias & ~choose_toward) = -1;

            % Unbiased +/-1 reorientation in the trap crop and BOTH
            % emigration zones. Each random draw maps to one actual index.
            idx = find(in_trap | outside);
            if ~isempty(idx)
                r = rand(numel(idx),1);
                dir_new(idx) = 1;
                dir_new(idx(r > 0.5)) = -1;
            end
            dir_new(~active) = 0;

            % Update positions of all active individuals.
            x(active) = x(active) + dir_new(active).*step_mag(active);
            dir = dir_new;

            % Entry interception requires crossing the interface FROM
            % THE CASH CROP. The endpoint may be in trap or beyond x=L.
            entered_now = was_in_cash & (x > zone_left);
            kill_on_entry = entered_now & ...
                            (rand(Npart,1) < gamma_entry);
            is_killed(kill_on_entry) = true;

            % Within-trap removal only applies to surviving individuals
            % whose new positions actually lie in the trap-crop strip.
            still_active = ~is_killed;
            in_trap_after_entry = still_active & ...
                                  (x > zone_left) & (x <= L);
            kill_inside = in_trap_after_entry & ...
                          (rand(Npart,1) < gamma_dt);
            is_killed(kill_inside) = true;

            % Permanently deactivate insects removed by either process.
            just_killed = kill_on_entry | kill_inside;
            x(just_killed) = NaN;
            dir(just_killed) = 0;
        end

        %% Record the four mutually exclusive states at day 15.
        alive = ~is_killed;
        nCash = sum(alive & (x >= 0) & (x <= zone_left));
        nTrap = sum(alive & (x > zone_left) & (x <= L));
        nRemoved = sum(is_killed);
        nEmigrated = sum(alive & ((x < 0) | (x > L)));

        assert(nCash+nTrap+nRemoved+nEmigrated == Npart, ...
            'Population accounting failed: scenario %d, k=%.2f.', ...
            icase, k_cash);

        cash_final(icase,ik) = 100*nCash/Npart;
        trap_final(icase,ik) = 100*nTrap/Npart;
        kill_final(icase,ik) = 100*nRemoved/Npart;
        emig_final(icase,ik) = 100*nEmigrated/Npart;
    end

    fprintf('Completed scenario %d of %d: %s\n', ...
        icase,nCases,labels{icase});
end

%% 5. VERIFY POPULATION CONSERVATION ACROSS ALL 36 SIMULATIONS
population_totals = cash_final+trap_final+kill_final+emig_final;
max_imbalance = max(abs(population_totals(:)-100));
assert(max_imbalance < 1e-10, ...
    'Four population states must sum to 100%% at every (scenario,k).');
fprintf('Maximum population-accounting discrepancy: %.3g percentage points\n', ...
    max_imbalance);

%% 6. PLOT FOUR PANELS (MATCHING FIGURE 6 COLOURS AND LINE STYLES)
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
hLegend = gobjects(nCases,1);

% (a) Percentage in cash crop at t=15 days.
axCash = subplot(2,2,1);
hold(axCash,'on');
for i = 1:nCases
    hLegend(i) = plot(k_vals,cash_final(i,:), ...
        'Color',colors(i,:),'LineStyle',lineStyles{i},'LineWidth',2);
end
ylim([0 100]);
xlim([k_vals(1) k_vals(end)]);
ylabel('Cash crop occupancy (%)');
title('(a) Cash crop occupancy','FontWeight','normal');
set(axCash,'FontName','Arial','FontSize',10, ...
    'Box','off','TickDir','out','LineWidth',1.5);
hold(axCash,'off');

% (b) Percentage in trap crop at t=15 days.
axTrap = subplot(2,2,2);
hold(axTrap,'on');
for i = 1:nCases
    plot(k_vals,trap_final(i,:), ...
        'Color',colors(i,:),'LineStyle',lineStyles{i},'LineWidth',2);
end
ylim([0 50]);
xlim([k_vals(1) k_vals(end)]);
ylabel('Trap crop occupancy (%)');
title('(b) Trap crop occupancy','FontWeight','normal');
set(axTrap,'FontName','Arial','FontSize',10, ...
    'Box','off','TickDir','out','LineWidth',1.5);
hold(axTrap,'off');

% (c) Cumulative removal at t=15 days.
axRemoved = subplot(2,2,3);
hold(axRemoved,'on');
for i = 1:nCases
    plot(k_vals,kill_final(i,:), ...
        'Color',colors(i,:),'LineStyle',lineStyles{i},'LineWidth',2);
end
ylim([0 100]);
xlim([k_vals(1) k_vals(end)]);
xlabel('Cash crop fraction, k','FontName','Arial');
ylabel('Removed (%)');
title('(c) Removed','FontWeight','normal');
set(axRemoved,'FontName','Arial','FontSize',10, ...
    'Box','off','TickDir','out','LineWidth',1.5);
hold(axRemoved,'off');

% (d) Current emigration at t=15 days (not cumulative departures).
axEmigrated = subplot(2,2,4);
hold(axEmigrated,'on');
for i = 1:nCases
    plot(k_vals,emig_final(i,:), ...
        'Color',colors(i,:),'LineStyle',lineStyles{i},'LineWidth',2);
end
ylim([0 50]);
xlim([k_vals(1) k_vals(end)]);
xlabel('Cash crop fraction, k','FontName','Arial');
ylabel('Emigrated (%)');
title('(d) Emigrated','FontWeight','normal');
set(axEmigrated,'FontName','Arial','FontSize',10, ...
    'Box','off','TickDir','out','LineWidth',1.5);
hold(axEmigrated,'off');

% Scenario names refer to the same colours/styles used in Figure 6.
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

fprintf('Figure 7 simulation complete.\n');
%% AUTOMATIC TIGHT PNG EXPORT (VISUAL OUTPUT ONLY)
% Removes blank canvas outside the plotted content, including labels,
% legend and colour bar. Does not change simulations or axes.
save_figure = true;    % Set false to display without exporting
if save_figure
    png_file = 'Fig_7.png';
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
