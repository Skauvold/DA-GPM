% A script for creating Fig. 9

clear all; close all; clc

% Load data
load('./data/thetaPriorRealisations.mat');
load('./data/thetapost3km');

% Set colors
priorcolor = [0.7, 0.7, 0.7];
postcolor = [0.3, 0.3, 0.3];

% Compute edges of confidence envelopes
slPriorLower = quantile(f_SL_ens', 0.05);
slPriorUpper = quantile(f_SL_ens', 0.95);

ssPriorLower = quantile(f_SS_ens', 0.05)/0.67;
ssPriorUpper = quantile(f_SS_ens', 0.95)/0.67;

slPostLower = quantile(theta_SL_final', 0.05);
slPostUpper = quantile(theta_SL_final', 0.95);

ssPostLower = quantile(theta_SS_final', 0.05)/0.67;
ssPostUpper = quantile(theta_SS_final', 0.95)/0.67;

% Create figure window and set subplot spacing parameters
figure()
gap = [0.02, 0.04];
marg_h = [0.12, 0.08];
marg_w = [0.12, 0.02];

% Set vertical axis limits
ylSS = [-0.3, 1.2];
ylSL = [-3000, 4000];

% Make upper left subplot
subtightplot(2,2,1, gap, marg_h, marg_w)
hold on
hss1 = fill([0:50, fliplr(0:50)], [ssPriorLower, fliplr(ssPriorUpper)], priorcolor);
hss2 = fill([0:50, fliplr(0:50)], [ssPostLower, fliplr(ssPostUpper)], postcolor);
hold off
alpha([hss1, hss2], 0.5)
set(gca, 'FontSize', 14)
ylabel('$\theta_\mathrm{SS}/\theta_\mathrm{SS,max}$', 'Interpreter', 'Latex', 'FontSize', 18)
[lh, lhicons] = legend([hss1,hss2], 'Prior', 'Posterior', 'Location', 'SouthWest');
set(lh, 'FontSize', 13);
PatchInLegend = findobj(lhicons, 'type', 'patch');
set(PatchInLegend, 'facealpha', 0.5)
box on
set(gca,'XTick',[])
xlim([0, 50])
ylim(ylSS)
title('$\sigma_{\Delta z} = 3$ km', 'Interpreter', 'Latex', 'FontSize', 18)

% Make lower left subplot
subtightplot(2,2,3, gap, marg_h, marg_w)
hold on
hsl1 = fill([0:50, fliplr(0:50)], [slPriorLower, fliplr(slPriorUpper)], priorcolor);
hsl2 = fill([0:50, fliplr(0:50)], [slPostLower, fliplr(slPostUpper)], postcolor);
hold off
alpha([hsl1, hsl2], 0.5)
xlabel('Time step', 'Interpreter', 'Latex', 'FontSize', 18);
set(gca, 'FontSize', 14)
ylabel('$\theta_\mathrm{SL}$ [m]', 'Interpreter', 'Latex', 'FontSize', 18)
[lh, lhicons] = legend([hsl1,hsl2], 'Prior', 'Posterior', 'Location', 'NorthEast');
set(lh, 'FontSize', 13);
PatchInLegend = findobj(lhicons, 'type', 'patch');
set(PatchInLegend, 'facealpha', 0.5)
box on
xlim([0, 50])
ylim(ylSL)

% Load more data
load('./data/thetapost30m');

% Compute more edges of confidence envelopes
slPostLower = quantile(theta_SL_final', 0.05);
slPostUpper = quantile(theta_SL_final', 0.95);

ssPostLower = quantile(theta_SS_final', 0.05)/0.67;
ssPostUpper = quantile(theta_SS_final', 0.95)/0.67;

% Make upper right plot
subtightplot(2,2,2, gap, marg_h, marg_w)
hold on
hss1 = fill([0:50, fliplr(0:50)], [ssPriorLower, fliplr(ssPriorUpper)], priorcolor);
hss2 = fill([0:50, fliplr(0:50)], [ssPostLower, fliplr(ssPostUpper)], postcolor);
hold off
alpha([hss1, hss2], 0.5)
set(gca, 'FontSize', 14)
[lh, lhicons] = legend([hss1,hss2], 'Prior', 'Posterior', 'Location', 'SouthWest');
set(lh, 'FontSize', 13);
PatchInLegend = findobj(lhicons, 'type', 'patch');
set(PatchInLegend, 'facealpha', 0.5)
box on
set(gca,'XTick',[], 'YTick',[])
xlim([0, 50])
ylim(ylSS)
title('$\sigma_{\Delta z} = 30$ m', 'Interpreter', 'Latex', 'FontSize', 18)

% Make lower right plot
subtightplot(2,2,4, gap, marg_h, marg_w)
hold on
hsl1 = fill([0:50, fliplr(0:50)], [slPriorLower, fliplr(slPriorUpper)], priorcolor);
hsl2 = fill([0:50, fliplr(0:50)], [slPostLower, fliplr(slPostUpper)], postcolor);
hold off
alpha([hsl1, hsl2], 0.5)
xlabel('Time step', 'Interpreter', 'Latex', 'FontSize', 18);
set(gca, 'FontSize', 14)
[lh, lhicons] = legend([hsl1,hsl2], 'Prior', 'Posterior', 'Location', 'NorthEast');
set(lh, 'FontSize', 13);
PatchInLegend = findobj(lhicons, 'type', 'patch');
set(PatchInLegend, 'facealpha', 0.5)
box on
set(gca, 'YTick',[])
xlim([0, 50])
ylim(ylSL)

% Set figure size and position
set(gcf, 'Units', 'Centimeters', 'Position', [10,10,19,13])