%% FIGURE 8 | ENTRY REMOVAL, WITHIN-TRAP REMOVAL AND ARRESTMENT
% Insect movement and spatial dynamics determine trap crop effectiveness
% Ecological Modelling (2026)
%
% PURPOSE
%   At t = 15 days, map cumulative removal (top row) and
%   CURRENT trap-crop occupancy (bottom row) against entry-triggered
%   removal gamma_entry and daily within-trap removal gamma_day.
%   Compare strong arrestment (rho = 0.01), baseline arrestment
%   (rho = 0.1), and no arrestment (rho = 1).
%
% MODEL GEOMETRY AND THE FOUR MUTUALLY EXCLUSIVE POPULATION STATES
%   Cash crop:        0 <= x <= k*L
%   Trap crop:        k*L < x <= L
%   Emigration:       x < 0 or x > L, with possible re-entry
%   Removed:          permanently inactive insects
%   All insects begin uniformly distributed within the cash crop.
%
% MOVEMENT RULES (IDENTICAL TO CORRECTED FIGURES 2-7)
%   At each step, x_i(t) = x_i(t-1) + delta_i(t)*s_i(t), where
%   s_i = sigma_region*abs(Z), Z ~ N(0,1).
%   In the cash crop, the direction mechanism is selected at each step:
%     Persistence with probability (1-omega):
%       P(same as previous direction) = (1+alpha)/2;
%     Attraction with probability omega:
%       beta(x) = beta_max*(x/(k*L))^nu,
%       P(rightward) = (1+beta(x))/2.
%   In the trap crop AND either emigration region, direction is a fresh,
%   unbiased +/-1 choice on each step. Only the trap-crop step-length
%   scale is reduced to sigma_trap = rho*sigma_cash; rho = 1 corresponds
%   to no arrestment. The emigrated insects remain active.
%
% REMOVAL RULES (IDENTICAL TO CORRECTED FIGURES 2-7)
%   On crossing FROM the cash crop past x = k*L, an individual is removed
%   with probability gamma_entry. This also includes the extremely rare
%   single step that overshoots the entire trap strip and ends at x > L.
%   There is no entry-removal trial upon returning to the trap strip from
%   right-hand emigration (x > L).
%   A surviving insect whose new position is INSIDE the trap crop is
%   removed with probability gamma_dt per step, where
%       gamma_dt = 1-(1-gamma_day)^(dt/24), dt in hours.
%   Thus gamma_day = 0 removes nobody via residence, and gamma_day = 1
%   removes every survivor present in the trap crop after a time step.
%
% SIMULATION DESIGN
%   One stochastic simulation of N = 5000 insects for each (rho,
%   gamma_day,gamma_entry) combination; results are not averaged over
%   stochastic replicates. A fixed seed permits reproducible figures.
%   Note that separate Figure 2-7 scripts consume random numbers in a
%   different order, so their stochastic outcomes need not match exactly
%   for the same parameter values.
%
% FIGURE LAYOUT
%   A 2x3 grid: top row cumulative removal, bottom row trap occupancy;
%   columns rho = 0.01, 0.1, and 1. All six panels use ONE shared
%   0-50%% colour scale. All six panels include labelled contours
%   at 10%%, 20%%, 30%%, and 40%%.
%
% COMPATIBILITY
%   MATLAB R2016a or later; no additional toolboxes required.

clear; clc; close all;
rng(1,'twister');

%% 1. FIXED MODEL PARAMETERS
N          = 5000;   % Initial insect population
L          = 100;    % Total field length (m)
k          = 0.8;    % Cash-crop fraction
zone_left  = k*L;    % Cash-crop/trap-crop boundary (m)
t_days     = 15;     % Snapshot time (days)
dt         = 0.5;    % Movement step duration (hours)
nSteps     = round(t_days*24/dt)+1; % Include the initial condition
sigma_cash = 1;      % Half-normal scale in cash and outside (m)

alpha      = 0.2;    % Persistence parameter in cash crop
omega      = 0.5;    % Chance of bias-driven movement in cash crop
beta_max   = 0.1;    % Maximum attraction strength at the trap boundary
nu         = 1;      % Linear attraction profile

%% 2. PARAMETER SWEEP
% gamma_entry acts at cash-to-trap entry; gamma_day acts during residence.
gamma_entry_vals = 0:0.01:0.25;
gamma_day_vals   = 0:0.01:0.25;
rho_vals         = [0.01,0.1,1];

nGE  = numel(gamma_entry_vals);
nGD  = numel(gamma_day_vals);
nRho = numel(rho_vals);

% Parameter checks: ensure all probabilities and movement scales are valid.
assert(N > 0 && L > 0 && k > 0 && k < 1 && dt > 0 && ...
       sigma_cash > 0,'Invalid population, geometry or step length.');
assert(abs(alpha) <= 1 && omega >= 0 && omega <= 1 && ...
       beta_max >= 0 && beta_max <= 1 && nu >= 0, ...
       'Invalid persistence or attraction parameters.');
assert(all(gamma_entry_vals >= 0 & gamma_entry_vals <= 1) && ...
       all(gamma_day_vals >= 0 & gamma_day_vals <= 1), ...
       'Removal probabilities must be between 0 and 1.');
assert(all(rho_vals > 0 & rho_vals <= 1), ...
       'Arrestment parameter must satisfy 0 < rho <= 1.');

%% 3. STORE FOUR POPULATION STATES AT 15 DAYS
% Dimensions: rho index, daily-removal index, entry-removal index.
final_cash_pct      = zeros(nRho,nGD,nGE);
final_trap_pct      = zeros(nRho,nGD,nGE);
final_removed_pct   = zeros(nRho,nGD,nGE);
final_emigrated_pct = zeros(nRho,nGD,nGE);

fprintf('Figure 8 scan: %d rho values x %d gamma_day values x %d gamma_entry values = %d runs.\n', ...
    nRho,nGD,nGE,nRho*nGD*nGE);

%% 4. SIMULATE EACH PARAMETER COMBINATION
for iR = 1:nRho
    rho = rho_vals(iR);
    sigma_trap = rho*sigma_cash;
    fprintf('\n=== rho = %.2f ===\n',rho);

    for iGD = 1:nGD
        gamma_day = gamma_day_vals(iGD);
        gamma_dt = 1-(1-gamma_day)^(dt/24);

        for iGE = 1:nGE
            gamma_entry = gamma_entry_vals(iGE);

            % Random initial positions and headings; all insects start
            % within the cash crop and remain active until removed.
            x = zone_left*rand(N,1);
            dir = sign(randn(N,1));
            dir(dir == 0) = 1;
            is_killed = false(N,1);

            for a = 2:nSteps
                active = ~is_killed;

                % CLASSIFY THE REGION AT THE START OF EACH STEP.
                in_cash = active & (x >= 0) & (x <= zone_left);
                in_trap = active & (x > zone_left) & (x <= L);
                outside = active & ((x < 0) | (x > L));
                was_in_cash = in_cash;

                % HALF-NORMAL STEP LENGTH; arrestment only within trap.
                sigma_eff = sigma_cash*ones(N,1);
                sigma_eff(in_trap) = sigma_trap;
                step_mag = abs(randn(N,1)).*sigma_eff;

                % SPATIALLY VARYING ATTRACTION, APPLIED ONLY IN CASH.
                beta_eff = zeros(N,1);
                beta_eff(in_cash) = beta_max*(x(in_cash)./zone_left).^nu;

                % PERSISTENCE-DRIVEN DIRECTION CHOICE IN CASH.
                dir_new = dir;
                keep_same = rand(N,1) < (1+alpha)/2;
                flip_idx = in_cash & ~keep_same;
                dir_new(flip_idx) = -dir_new(flip_idx);

                % INDEPENDENTLY CHOOSE BIAS-DRIVEN MOVEMENT IN CASH.
                do_bias = in_cash & (rand(N,1) < omega);
                choose_toward = rand(N,1) < (1+beta_eff)/2;
                dir_new(do_bias) = 1;
                dir_new(do_bias & ~choose_toward) = -1;

                % UNBIASED SIMPLE RANDOM WALK IN TRAP AND BOTH EMIGRATION
                % REGIONS. FIX: FIND EXPLICIT INDICES BEFORE SUBINDEXING.
                idx = find(in_trap | outside);
                if ~isempty(idx)
                    r = rand(numel(idx),1);
                    dir_new(idx) = 1;
                    dir_new(idx(r > 0.5)) = -1;
                end
                dir_new(~active) = 0;

                % UPDATE POSITIONS OF ACTIVE INSECTS ONLY.
                x(active) = x(active) + dir_new(active).*step_mag(active);
                dir = dir_new;

                % ENTRY REMOVAL: CASH -> TRAP INTERFACE, INCLUDING THE
                % EXTREMELY RARE OVERSHOOT FROM CASH TO BEYOND x=L.
                entered_now = was_in_cash & (x > zone_left);
                kill_on_entry = entered_now & ...
                                (rand(N,1) < gamma_entry);
                is_killed(kill_on_entry) = true;

                % WITHIN-TRAP REMOVAL AMONG SURVIVING TRAP RESIDENTS.
                still_active = ~is_killed;
                in_trap_after_entry = still_active & ...
                                      (x > zone_left) & (x <= L);
                kill_inside = in_trap_after_entry & ...
                              (rand(N,1) < gamma_dt);
                is_killed(kill_inside) = true;

                % Removed insects never take another step.
                just_killed = kill_on_entry | kill_inside;
                x(just_killed) = NaN;
                dir(just_killed) = 0;
            end

            %% Count the four mutually exclusive states.
            alive = ~is_killed;
            nCash = sum(alive & (x >= 0) & (x <= zone_left));
            nTrap = sum(alive & (x > zone_left) & (x <= L));
            nRemoved = sum(is_killed);
            nEmigrated = sum(alive & ((x < 0) | (x > L)));

            assert(nCash+nTrap+nRemoved+nEmigrated == N, ...
                'Population accounting failed: rho=%g, gamma_day=%g, gamma_entry=%g.', ...
                rho,gamma_day,gamma_entry);

            % Exact limiting-case checks for both removal mechanisms.
            % With gamma_day=1, all insects ending a step inside the
            % trap are killed. Perfect entry interception also leaves
            % no living insects in the trap when all insects start in cash.
            if gamma_day == 1 || gamma_entry == 1
                assert(nTrap == 0, ...
                    'Trap occupancy must be zero under perfect removal.');
            end
            if gamma_entry == 0 && gamma_day == 0
                assert(nRemoved == 0, ...
                    'No removal is possible when both mechanisms are off.');
            end

            final_cash_pct(iR,iGD,iGE) = 100*nCash/N;
            final_trap_pct(iR,iGD,iGE) = 100*nTrap/N;
            final_removed_pct(iR,iGD,iGE) = 100*nRemoved/N;
            final_emigrated_pct(iR,iGD,iGE) = 100*nEmigrated/N;
        end

        if mod(iGD,2) == 0
            fprintf('  Finished gamma_day = %.2f\n', gamma_day);
        end
    end
end

%% 5. CHECK POPULATION CONSERVATION FOR EVERY PARAMETER COMBINATION
total_pct = final_cash_pct + final_trap_pct + ...
            final_removed_pct + final_emigrated_pct;
max_imbalance = max(abs(total_pct(:)-100));
assert(max_imbalance < 1e-10, ...
    'Four population-state percentages must sum to 100%%.');
fprintf('\nMaximum population-accounting discrepancy: %.3g percentage points\n', ...
    max_imbalance);

% Report the data ranges for checking the shared 0-50% colour scale.
fprintf('Removed range:  %.2f to %.2f%%\n', ...
    min(final_removed_pct(:)),max(final_removed_pct(:)));
fprintf('Trap occupancy range: %.2f to %.2f%%\n', ...
    min(final_trap_pct(:)),max(final_trap_pct(:)));
fprintf('Emigrated range (not plotted): %.2f to %.2f%%\n', ...
    min(final_emigrated_pct(:)),max(final_emigrated_pct(:)));

%% 6. ONE SHARED COLOUR SCALE FOR ALL SIX PANELS
% Fixed at 0-50% as requested; no model outputs are changed.
common_lim = [0 50];
max_plotted = max([final_removed_pct(:); final_trap_pct(:)]);
if max_plotted > common_lim(2)
    warning('Values reach %.2f%%; the shared 0-50%% colour scale clips values above 50%%.', max_plotted);
end
fprintf('Shared colour range: %.1f to %.1f%%\n',common_lim);

%% 7. PUBLICATION-STYLE THREE-COLUMN / TWO-ROW FIGURE
% Columns: rho = 0.01, 0.1, 1; rows: removal, trap occupancy.
figure('Color','w','Position',[50 50 1450 840]);
colormap(jet(256));

% Explicit positions leave room for the single shared colour bar.
x_positions = [0.065 0.355 0.645];
y_positions = [0.575 0.155];  % Removal, then trap occupancy
ax_width = 0.220;
ax_height = 0.330;
ax_positions = zeros(6,4);
ax = gobjects(6,1);
panel_labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};

for iOutcome = 1:2
    for iR = 1:nRho
        panel = (iOutcome-1)*nRho + iR;
        ax_positions(panel,:) = [x_positions(iR),y_positions(iOutcome),ax_width,ax_height];
        ax(panel) = axes('Position',ax_positions(panel,:));

        if iOutcome == 1
            Z = squeeze(final_removed_pct(iR,:,:));
            out_title = 'Removed';
        else
            Z = squeeze(final_trap_pct(iR,:,:));
            out_title = 'Trap-crop occupancy';
        end

        imagesc(gamma_entry_vals,gamma_day_vals,Z);
        set(ax(panel),'YDir','normal');
        caxis(ax(panel),common_lim);
        hold(ax(panel),'on');

        % Use the same contour levels in all panels for consistency.
        contour_levels = 10:10:(common_lim(2)-0.01);
        if ~isempty(contour_levels)
            [C,h] = contour(ax(panel),gamma_entry_vals, ...
                gamma_day_vals,Z,contour_levels, ...
                'LineColor','k','LineWidth',0.7);
            clabel(C,h,'FontName','Arial','FontSize',8, ...
                'Color','k','LabelSpacing',400);
        end
        hold(ax(panel),'off');

        set(ax(panel),'FontName','Arial','FontSize',10, ...
            'Box','on','LineWidth',1,'TickDir','out', ...
            'XLim',[0 0.25],'YLim',[0 0.25]);
        axis(ax(panel),'square');
        title(ax(panel),sprintf('%s %s, \\rho = %g', ...
            panel_labels{panel},out_title,rho_vals(iR)), ...
            'FontName','Arial','FontSize',12, ...
            'FontWeight','normal','Interpreter','tex');

        if iOutcome == 2
            xlabel(ax(panel), ...
                '\gamma_{entry} (entry-triggered removal probability)', ...
                'FontName','Arial','FontSize',10,'Interpreter','tex');
        end
        if iR == 1
            ylabel(ax(panel), ...
                '\gamma_{day} (within-trap crop removal probability)', ...
                'FontName','Arial','FontSize',10,'Interpreter','tex');
        end
    end
end

%% 8. SINGLE SHARED COLOUR BAR FOR ALL SIX PANELS
cb = colorbar('peer',ax(6));
set(cb,'Position',[0.920 0.155 0.016 0.735], ...
    'FontName','Arial','FontSize',10);
cb.Label.String = 'Proportion of initial population (%)';
cb.Label.FontName = 'Arial';
cb.Label.FontSize = 10;
cb.Ticks = 0:10:50;

% Restore subplot positions after attaching the colour bar.
for p = 1:6
    set(ax(p),'Position',ax_positions(p,:));
end

fprintf('\nFigure 8 complete: 2x3 layout, rho = 0.01, 0.1, 1.\n');
%% AUTOMATIC TIGHT PNG EXPORT (VISUAL OUTPUT ONLY)
% Removes blank canvas outside the plotted content, including labels,
% legend and colour bar. Does not change simulations or axes.
save_figure = true;    % Set false to display without exporting
if save_figure
    png_file = 'Fig_8.png';
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
