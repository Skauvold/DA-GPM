% A script for creating Fig. 4

clear all; close all; clc

% Time axis
tt = 0:20;

sel = 3;

% Set vertical axis limits
yl1 = [-5, 35];
yl2 = [-40, 80];

% Set subplot parameters
hmarg = [0.12, 0.08];
wmarg = [0.10, 0.04];
gap = [0.08, 0.05];

% Create figure window
figure();

% Load data
load('./data/LUTSSLL_EnKF.mat')

% Create upper left plot
subtightplot(2,3,1, gap, hmarg, wmarg)
hold on
filh = fill([tt,fliplr(tt)],[LSS;flipud(USS)],0.8*[1,1,1]);
set(filh,'edgecolor','none')
plot(tt,TSS,'k-','LineWidth',1);
title('EnKF');
hold off
box on
ylabel('$\theta_\mathrm{SS}$ [mm/yr]','Interpreter','Latex','FontSize',18);
set(gca,'XtickLabel',[],'FontSize',14, 'YTick', [0, 40/4, 80/4, 120/4], 'YTickLabel', [0, 40, 80, 120]);
ylim(yl1)

% Create lower left plot
subtightplot(2,3,4, gap, hmarg, wmarg)
hold on
filh = fill([tt,fliplr(tt)],[LSL;flipud(USL)],0.8*[1,1,1]);
set(filh,'edgecolor','none')
plot(tt,TSL,'k-','LineWidth',1);
hold off
box on
xlabel('Time step', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$\theta_\mathrm{SL}$ [m]','Interpreter','Latex','FontSize',18);
set(gca,'FontSize',14, 'YTick', [-40, 0, 40, 80])
ylim(yl2)

% Create upper middle plot
load('./data/LUTSSLL_EnS.mat')
subtightplot(2,3,2, gap, hmarg, wmarg)
hold on
filh = fill([tt,fliplr(tt)],[LSS;flipud(USS)],0.8*[1,1,1]);
set(filh,'edgecolor','none')
plot(tt,TSS,'k-','LineWidth',1);
title('EnS');
hold off
box on
set(gca,'XtickLabel',[],'FontSize',14);
set(gca,'YTickLabel',[])
ylim(yl1)

% Create lower middle plot
subtightplot(2,3,5, gap, hmarg, wmarg)
hold on
filh = fill([tt,fliplr(tt)],[LSL;flipud(USL)],0.8*[1,1,1]);
set(filh,'edgecolor','none')
plot(tt,TSL,'k-','LineWidth',1);
hold off
box on
xlabel('Time step', 'Interpreter', 'Latex', 'FontSize', 18)
set(gca,'FontSize',14);
set(gca,'YTickLabel',[])
ylim(yl2)

% Load more data
load('./data/LUTSSLL_MDA.mat')

% Compute confidence envelope edges
LSS = quantile(SSres', 0.1)';
USS = quantile(SSres', 0.9)';
TSS = SSref;

LSL = quantile(SLres', 0.1)';
USL = quantile(SLres', 0.9)';
TSL = SLref;

% Create upper right plot
subtightplot(2,3,3, gap, hmarg, wmarg)
hold on
filh = fill([tt,fliplr(tt)],[LSS;flipud(USS)],0.8*[1,1,1]);
set(filh,'edgecolor','none')
plot(tt,TSS,'k-','LineWidth',1);
title('MDA');
hold off
box on
set(gca,'XtickLabel',[],'FontSize',14);
set(gca,'YTickLabel',[])
ylim(yl1)

% Create lower right plot
subtightplot(2,3,6, gap, hmarg, wmarg)
hold on
filh = fill([tt,fliplr(tt)],[LSL;flipud(USL)],0.8*[1,1,1]);
set(filh,'edgecolor','none')
plot(tt,TSL,'k-','LineWidth',1);
hold off
box on
xlabel('Time step', 'Interpreter', 'Latex', 'FontSize', 18)
set(gca,'FontSize',14);
set(gca,'YTickLabel',[])
ylim(yl2)

set(gcf, 'units', 'centimeters', 'position', [15, 30, 18, 12]);