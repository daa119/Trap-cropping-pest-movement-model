%% FIGURE 4 | CASH-CROP OCCUPANCY ACROSS MOVEMENT PERSISTENCE AND BIAS
% Insect movement and spatial dynamics determine trap crop effectiveness
% Ecological Modelling (2026)
%
% PURPOSE
%   Evaluate the proportion of insects remaining in the cash crop after
%   15 days as a function of the persistence parameter alpha and the
%   frequency omega of selecting environmentally biased movement.
%   Four panels correspond to beta_max = 0, 0.1, 0.25 and 0.5.
%   Each panel displays omega = 0, 0.25, 0.5, 0.75 and 1.
%
% MODEL GEOMETRY
%   The one-dimensional field occupies [0,L]. Its cash crop occupies
%   [0,k*L], and its terminal trap crop occupies (k*L,L]. Insects can
%   temporarily emigrate into x<0 or x>L and later re-enter the field.
%
% MOVEMENT PROCESS
%   A surviving insect's position is updated by
%       x_i(t) = x_i(t-1) + delta_i(t)*s_i(t),
%   with half-normal step length s_i = sigma_eff*abs(Z), Z~N(0,1).
%   In the cash crop, direction is selected using
%       P(use attraction) = omega;
%       P(use persistence) = 1-omega.
%   The persistence rule maintains the previous direction with probability
%       P(same heading) = (1+alpha)/2.
%   Attraction favours movement towards the trap crop (positive x):
%       beta(x) = beta_max*(x/(k*L))^nu;
%       P(rightward | attraction) = (1+beta(x))/2.
%   Movement is unbiased inside the trap crop and outside the field:
%       P(leftward) = P(rightward) = 1/2.
%   The trap-crop step-length scale is reduced to rho*sigma_cash.
%
% REMOVAL PROCESS
%   Insects entering the trap crop from the cash crop are removed with
%   probability gamma_entry, including on subsequent re-entries.
%   Insects that remain in the trap crop are independently subject to
%   per-step removal with probability
%       gamma_dt = 1 - (1-gamma_day)^(dt/24), for dt measured in hours.
%
% SIMULATION DESIGN AND REPRODUCIBILITY
%   One N=5000-particle stochastic realisation is simulated for each
%   (beta_max,omega,alpha) parameter combination. The complete sweep is
%   used to evaluate the desired combinations and plot the five selected
%   omega values. There is no averaging over simulation replicates.
%   The random-number seed is fixed for reproducibility.
%
% COMPATIBILITY
%   MATLAB R2016a or later; no additional toolboxes required.

clear; clc; close all;
rng(1, 'twister');

%% 1. FIXED PARAMETERS (BASELINE FROM FIGURE 2)
N           = 5000;                     % Initial number of insects
L           = 100;                      % Field length (m)
k           = 0.8;                      % Cash-crop fraction
t_days      = 15;                       % Evaluation time (days)
dt          = 0.5;                      % Movement step duration (h)
nSteps      = round(t_days*24/dt) + 1;  % Includes initial time point
zone_left   = k*L;                      % Cash/trap interface (m)

sigma       = 1;                        % Cash-crop half-normal scale (m)
rho         = 0.1;                      % Arrestment factor
sigma_trap  = rho*sigma;                % Trap-crop half-normal scale (m)

gamma_entry = 0.1;                      % Removal probability on entry
gamma_day   = 0.2;                      % Removal probability over 24 h
gamma_dt    = 1 - (1-gamma_day)^(dt/24);% Per-step removal probability
nu          = 1;                        % Linear spatial attraction profile

%% 2. MOVEMENT-PARAMETER SWEEP
Bmax_vals  = [0 0.1 0.25 0.5];          % Four panels
alpha_vals = -1:0.1:1;                  % Correlation between successive steps
omega_vals = 0:0.05:1;                  % Frequency of cue-driven movement

nBmax = numel(Bmax_vals);
nA    = numel(alpha_vals);
nO    = numel(omega_vals);
cash_occ_ALL = nan(nA,nO,nBmax);

% Parameter-domain checks.
assert(N > 0 && L > 0 && dt > 0 && sigma > 0,'Positive model parameters required.');
assert(k > 0 && k < 1,'Cash-crop fraction must satisfy 0 < k < 1.');
assert(rho > 0 && rho <= 1,'Require 0 < rho <= 1.');
assert(all(Bmax_vals >= 0 & Bmax_vals <= 1),'Require 0 <= beta_max <= 1.');
assert(all(abs(alpha_vals) <= 1),'Require -1 <= alpha <= 1.');
assert(all(omega_vals >= 0 & omega_vals <= 1),'Require 0 <= omega <= 1.');
assert(nu >= 0,'Require nu >= 0.');
assert(gamma_entry >= 0 && gamma_entry <= 1 && ...
       gamma_day >= 0 && gamma_day <= 1,'Removal probabilities must lie within [0,1].');

fprintf('Running alpha-omega scan for %d beta_max values ...\n',nBmax);
max_population_imbalance = 0;

%% 3. STOCHASTIC SIMULATIONS
for iBmax = 1:nBmax
    beta_max = Bmax_vals(iBmax);
    fprintf('\n=== beta_max = %.2f ===\n',beta_max);

    for jO = 1:nO
        omega = omega_vals(jO);

        for iA = 1:nA
            alpha = alpha_vals(iA);

            % Uniform initial positions and independent initial headings.
            x = zone_left*rand(N,1);
            dir = sign(randn(N,1));
            dir(dir==0) = 1;
            is_killed = false(N,1);

            for a = 2:nSteps
                active = ~is_killed;

                % Determine each insect's region before the step.
                in_cash = active & (x >= 0) & (x <= zone_left);
                in_trap = active & (x > zone_left) & (x <= L);
                outside = active & ((x < 0) | (x > L));
                was_in_cash = in_cash;

                % Half-normal step lengths, with arrestment only in trap.
                sigma_eff = sigma*ones(N,1);
                sigma_eff(in_trap) = sigma_trap;
                step_mag = abs(randn(N,1)).*sigma_eff;

                % Attraction strength increases towards the interface.
                beta_eff = zeros(N,1);
                beta_eff(in_cash) = beta_max*...
                    (x(in_cash)./zone_left).^nu;

                % Persistence-driven directions in the cash crop.
                dir_new = dir;
                keep_same = rand(N,1) < (1+alpha)/2;
                reverse = in_cash & ~keep_same;
                dir_new(reverse) = -dir_new(reverse);

                % Independently decide which insects use environmental bias.
                use_bias = in_cash & (rand(N,1) < omega);
                choose_right = rand(N,1) < (1+beta_eff)/2;
                dir_new(use_bias) = 1;
                dir_new(use_bias & ~choose_right) = -1;

                % Unbiased direction choices inside trap and outside field.
                idx = find(in_trap | outside);
                if ~isempty(idx)
                    r = rand(numel(idx),1);
                    dir_new(idx) = 1;
                    dir_new(idx(r > 0.5)) = -1;
                end
                dir_new(~active) = 0;

                % Position update.
                x(active) = x(active) + dir_new(active).*step_mag(active);
                dir = dir_new;

                % Entry-triggered removal when crossing the cash/trap interface.
                % This also catches a rare step crossing the entire strip.
                entered_now = was_in_cash & (x > zone_left);
                kill_entry = entered_now & (rand(N,1) < gamma_entry);
                is_killed(kill_entry) = true;
                x(kill_entry) = NaN;
                dir(kill_entry) = 0;

                % Within-trap removal among surviving trap-crop occupants.
                still_active = ~is_killed;
                in_trap_now = still_active & (x > zone_left) & (x <= L);
                kill_trap = in_trap_now & (rand(N,1) < gamma_dt);
                is_killed(kill_trap) = true;
                x(kill_trap) = NaN;
                dir(kill_trap) = 0;
            end

            % Evaluate cash-crop occupancy after 15 days.
            alive = ~is_killed;
            nCash = sum(alive & (x >= 0) & (x <= zone_left));
            nTrap = sum(alive & (x > zone_left) & (x <= L));
            nRemoved = sum(is_killed);
            nEmigrated = sum(alive & ((x < 0) | (x > L)));

            % Diagnostic: every insect occupies exactly one population state.
            imbalance = abs(nCash+nTrap+nRemoved+nEmigrated-N);
            max_population_imbalance = max(max_population_imbalance,imbalance);
            assert(imbalance == 0,'Population conservation failed.');

            cash_occ_ALL(iA,jO,iBmax) = 100*nCash/N;
        end

        if mod(jO,5)==0
            fprintf('  Finished omega = %.2f\n',omega);
        end
    end
end

fprintf('\nSimulation complete: %d parameter combinations.\n',nBmax*nO*nA);
fprintf('Maximum population-count imbalance: %g insects.\n',max_population_imbalance);

%% 4. FIGURE 4 — FOUR PANELS OF CASH-CROP OCCUPANCY
% The full omega sweep uses increments of 0.05. Display only the five
% selected omega curves defined in the manuscript's figure caption.
omega_plot_vals = [0 0.25 0.5 0.75 1];

% Consistent 2-by-2 plot-box dimensions and spacing, matched to Figure 2.
% Units are pixels; all four axes have the same width and height.
% This is strictly a plotting adjustment; the simulation is unchanged.
panel_positions_px = [ ...
     70 460 335 290; ...  % (a) upper left
    535 460 335 290; ...  % (b) upper right
     70  65 335 290; ...  % (c) lower left
    535  65 335 290];     % (d) lower right

figure('Color','w','Position',[60 60 900 800]);
% Four MATLAB line patterns plus a distinct marker-only fifth style.
% MATLAB does not provide a fifth connected LineStyle.
lineStyles = {'-','--',':','-.','none'};
markers    = {'o','s','^','d','x'};
markerSize = 2.2;
lineWidth  = 2;

panel_labels = {'(a)','(b)','(c)','(d)'};
panel_axes = gobjects(4,1);

for iBmax = 1:nBmax
    ax = subplot(2,2,iBmax);
    panel_axes(iBmax) = ax;
    hold(ax,'on');

    cash_occ_AO = cash_occ_ALL(:,:,iBmax);
    for jP = 1:numel(omega_plot_vals)
        % Select the simulated omega value closest to the desired curve.
        [~,idxO] = min(abs(omega_vals-omega_plot_vals(jP)));
        ls = lineStyles{jP};
        mk = markers{jP};

        if strcmp(ls,'none')
            % Fifth omega curve: distinct x markers at every simulated alpha.
            % All plotted values come directly from cash_occ_AO, unchanged.
            plot(alpha_vals,cash_occ_AO(:,idxO), ...
                'LineStyle','none','Marker',mk,'MarkerSize',7, ...
                'LineWidth',1.5,'Color','k', ...
                'DisplayName',sprintf('\\omega = %.2f',omega_vals(idxO)));
        else
            plot(alpha_vals,cash_occ_AO(:,idxO), ...
                'LineStyle',ls,'Color','k','LineWidth',lineWidth, ...
                'DisplayName',sprintf('\\omega = %.2f',omega_vals(idxO)));

            % Markers are interpolated along the displayed polyline for clarity.
            alpha_marker_vals = -1:0.25:1;
            cash_marker_vals = interp1(alpha_vals,cash_occ_AO(:,idxO), ...
                alpha_marker_vals,'linear');
            plot(alpha_marker_vals,cash_marker_vals, ...
                'LineStyle','none','Marker',mk,'MarkerSize',markerSize, ...
                'MarkerFaceColor','w','MarkerEdgeColor','k', ...
                'HandleVisibility','off');
        end
    end

    if iBmax >= 3
        xlabel('\alpha (directional persistence)', ...
            'FontSize',11,'FontName','Arial');
    end
    if iBmax==1 || iBmax==3
        ylabel('Cash-crop occupancy (%)', ...
            'FontSize',11,'FontName','Arial');
    end

    title(sprintf('%s \\beta_{max} = %g',panel_labels{iBmax},Bmax_vals(iBmax)), ...
        'FontSize',12,'FontName','Arial','FontWeight','normal', ...
        'Interpreter','tex');
    grid off;
    ylim([0 100]);
    xlim([min(alpha_vals) max(alpha_vals)]);
    set(gca,'XTick',-1:0.25:1);
    set(gca,'FontName','Arial','FontSize',10,'Box','off', ...
        'TickDir','out','LineWidth',1.5);
    hold(ax,'off');
end

% A legend in the final panel differentiates the five omega settings.
ax4 = subplot(2,2,4);
hold(ax4,'on');
hLegendD = gobjects(numel(omega_plot_vals),1);
for jP = 1:numel(omega_plot_vals)
    if strcmp(lineStyles{jP},'none')
        % Display the fifth marker-only curve correctly in the legend.
        hLegendD(jP) = plot(ax4,nan,nan, ...
            'LineStyle','none','Marker',markers{jP}, ...
            'MarkerSize',7,'LineWidth',1.5,'Color','k', ...
            'DisplayName',sprintf('\\omega = %.2f',omega_plot_vals(jP)));
    else
        hLegendD(jP) = plot(ax4,nan,nan, ...
            'LineStyle',lineStyles{jP},'Color','k', ...
            'LineWidth',lineWidth, ...
            'DisplayName',sprintf('\\omega = %.2f',omega_plot_vals(jP)));
    end
end
legend(ax4,hLegendD,get(hLegendD,'DisplayName'), ...
    'Location','northeast','Box','off', ...
    'FontSize',10,'FontName','Arial');
hold(ax4,'off');

% Lock the plotting boxes to the common four-panel layout after legends
% (and colour bars, if present) have been created.
for iPanel = 1:4
    set(panel_axes(iPanel), 'Units', 'pixels', ...
        'Position', panel_positions_px(iPanel,:));
end
%% AUTOMATIC TIGHT PNG EXPORT (VISUAL OUTPUT ONLY)
% Removes blank canvas outside the plotted content, including labels,
% legend and colour bar. Does not change simulations or axes.
save_figure = true;    % Set false to display without exporting
if save_figure
    png_file = 'Fig_4.png';
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
