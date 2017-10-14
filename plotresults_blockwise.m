% Load, post-process and visualize results of North Slope, Alaska data
% assimilation case
%
% For Fig. 7, see l. 201.
% For Fig. 8, see l. 423.

clear; close all; clc

% Load results
load('./data/ANSCaseEnKF_output.mat');

%% Processing

% Time axis properties
steps = 50;
tlen = 51;

% Ensemble size
ne = 100;

final_analysis = Xa;

theta_SL_final = final_analysis(1:tlen,:);
theta_SS_final = final_analysis(tlen+1:2*tlen,:);
ZS_final = final_analysis(2*tlen+1:end,:);

% Size of model grid
ny = 87;
nx = 110;

%% Make GR conditioning figure

numberOfBlocks = 10; % 10 blocks of 5 layers each

% Step through blocks and rearrange data for plotting
for block = 1:numberOfBlocks
    grobs(block,:) = observed_data{block}(2,:);
    grapr(block,:) = synthetic_apriori{block}(2,:);
    grapo(block,:) = synthetic_aposteriori{block}(2,:);
    graprapo(2*block-1,:) = grapr(block,:);
    graprapo(2*block,:) = grapo(block,:);
end

% Define axis for plotting
blockaxis = 5:5:50;
blockaxis = [blockaxis;blockaxis];
blockaxis = [0;blockaxis(:);55]';
blockaxis2 = [2.5:5:47.5]';

% Make figures
figure()
hold on
for i = 1:10
    plot(5*[i,i], [30,100], 'k--');
end
plot(blockaxis, [graprapo(1,:); graprapo; graprapo(end,:)], '-', 'Color', 0.75*ones(1,3));
hold off
xlabel('Layer, k')
ylabel('GR [API]')
title('GR log fit, Blockwise EnKF')

figure()
hold on
stairs(0:numberOfBlocks, [grapr;grapr(end,:)], 'Color', 0.75*ones(1,3));
plot(0.5:(numberOfBlocks-0.5), grobs, 'k.');
plot(0.5:(numberOfBlocks-0.5), mean(grobs, 2),'k.', 'MarkerSize', 25);
hold off
xlabel('Block number')
ylabel('GR [API]')
title('A priori')

figure()
hold on
stairs(0:numberOfBlocks, [grapo;grapo(end,:)], 'Color', 0.75*ones(1,3));
plot(0.5:(numberOfBlocks-0.5), grobs, 'k.');
plot(0.5:(numberOfBlocks-0.5), mean(grobs, 2),'k.', 'MarkerSize', 25);
hold off
xlabel('Block number')
ylabel('GR [API]')
title('A posteriori')

%% Plot theta_SL and theta_SS

% Load parameter prior ensemble
load('./data/thetaPriorRealisations.mat');

figure()
hold on
plot(0:50, f_SL_ens, 'b-')
plot(0:50, theta_SL_final, 'r-')
hold off
title('Sea level, 100 ensemble members')
xlabel('Time step')
set(gca, 'FontSize', 14)
ylabel('$\theta_\mathrm{SL}$', 'Interpreter', 'Latex', 'FontSize', 18)
box on

figure()
hold on
plot(0:50, f_SS_ens, 'b-')
plot(0:50, theta_SS_final, 'r-')
hold off
title('Sediment supply, 100 ensemble members')
xlabel('Time step')
set(gca, 'FontSize', 14)
ylabel('$\theta_\mathrm{SS}$', 'Interpreter', 'Latex', 'FontSize', 18)
box on

% Define colours
lightgreen = [0.3, 1.0, 0.3];
lightblue = [0.0, 0.3, 0.4];
lightred = [1.0, 0.3, 0.3];

% Compute confidence envelope edges
slPriorLower = quantile(f_SL_ens', 0.05);
slPriorUpper = quantile(f_SL_ens', 0.95);

slPostLower = quantile(theta_SL_final', 0.05);
slPostUpper = quantile(theta_SL_final', 0.95);

ssPriorLower = quantile(f_SS_ens', 0.05)/0.67;
ssPriorUpper = quantile(f_SS_ens', 0.95)/0.67;

ssPostLower = quantile(theta_SS_final', 0.05)/0.67;
ssPostUpper = quantile(theta_SS_final', 0.95)/0.67;

figure()
hold on
hsl1 = fill([0:50, fliplr(0:50)], [slPriorLower, fliplr(slPriorUpper)], lightblue);
hsl2 = fill([0:50, fliplr(0:50)], [slPostLower, fliplr(slPostUpper)], lightred);
hold off
title('Sea level, 90% confidence')
alpha([hsl1, hsl2], 0.5)
xlabel('Time step')
set(gca, 'FontSize', 14)
ylabel('$\theta_\mathrm{SL}$ [m]', 'Interpreter', 'Latex', 'FontSize', 18)
legend([hsl1,hsl2], 'Prior', 'Posterior', 'Location', 'SouthEast')
box on

figure()
hold on
hss1 = fill([0:50, fliplr(0:50)], [ssPriorLower, fliplr(ssPriorUpper)], lightblue);
hss2 = fill([0:50, fliplr(0:50)], [ssPostLower, fliplr(ssPostUpper)], lightred);
hold off
title('Sediment supply, 90% confidence')
alpha([hss1, hss2], 0.5)
xlabel('Time step')
set(gca, 'FontSize', 14)
ylabel('$\theta_\mathrm{SS}/\theta_\mathrm{SS,max}$', 'Interpreter', 'Latex', 'FontSize', 18)
legend([hss1,hss2], 'Prior', 'Posterior', 'Location', 'SouthWest')
box on

%% Vertical plot of synthetic GR observations

figure()
hold on

grValues = [25, 40, 55, 115];

xw = 67-2;
yw = 69-3;

zwArray = [];
grArray = [];

for b = 1:ne
    ZS_array = reshape(ZS_final(:,b),nx,ny,[]);
    Z_final = ZS_array(:,:,1:steps+1);
    S1_final = ZS_array(:,:,steps+2:2*steps+1);
    S2_final = ZS_array(:,:,2*steps+2:3*steps+1);
    S3_final = ZS_array(:,:,3*steps+2:4*steps+1);
    
    zw = [];
    grw = [];
    
    zwArray(1,b) = Z_final(xw,yw,1);
    
    for i = 1:steps
        zw(i) = Z_final(xw,yw,i+1);
        zwArray(i+1,b) = zw(i);
        s1wi = S1_final(xw,yw,i);
        s2wi = S2_final(xw,yw,i);
        s3wi = S3_final(xw,yw,i);
        zgrwi = gr_obs_model([zw(i);s1wi;s2wi;s3wi],grValues);
        grw(i) = zgrwi(2);
        grArray(i,b) = grw(i);
    end
    plotHandle1 = stairs([grw(1),grw(1:end-1)],zw,'r-','LineWidth',1);
end

xlabel('Tunalik 1 GR, Updated','FontSize',14)
ylabel('z [m]','FontSize',14)
set(gca,'XAxislocation','Top','FontSize',14)
box on

load('./data/Tunalik1data.mat');

xlabel('$\gamma$','Interpreter','Latex','FontSIze',18)
ylabel('z')
set(gca,'FontSize',14)

%% Comparison of prior and posterior thickness (Fig. 7)

% Initialise arrays
zCombinedArray = [];
zPriorArray = [];
zPostArray = [];
zObservationsArray = [];

% Step through blocks and prepare data for plotting
for k = 1:10
    zPriorArray(k,:) = synthetic_apriori{k}(1,:);
    zPostArray(k,:) = synthetic_aposteriori{k}(1,:);
    zObservationsArray(k,:) = observed_data{k}(1,:);
    zCombinedArray(2*k-1,:) = synthetic_apriori{k}(1,:);
    zCombinedArray(2*k,:) = synthetic_aposteriori{k}(1,:);
end

zCombinedArray = [zeros(1,ne); zCombinedArray; zCombinedArray(end,:)];

% Define colour for ensemble members
lightgray = 0.67*ones(1,3);

% Define blockwise time axis for plotting
kk = repmat([5:5:50],2,1);
kk = kk(:);
kk = [0,kk'];

% Create figure window
figure()

% Set parameters for subplot spacing
gap = [0.05, 0.05];
marg_h = [0.15, 0.10];
marg_w = [0.15, 0.05];
yl = [0, 2400];

% Create lower left subplot
subtightplot(2,2,3, gap, marg_h, marg_w)
dzArrayPost = zwArray - zwArray(1,:);
load('./data/dzArrayPrior.mat');
dzArrayPrior = [zeros(1, size(dzArray,2)); dzArray];
hold on
priorHandle = plot(0:50, dzArrayPrior, '-', 'Color', 0.7*ones(1,3));
posteriorHandle = plot(0:50, dzArrayPost, 'k-');
hold off
box on
legendHandle = legend([priorHandle(1), posteriorHandle(1)], 'Prior', 'Posterior', 'Location', 'NorthWest');
set(gca,'XTick',[0, 10, 20, 30, 40, 50], 'FontSize', 14)
set(legendHandle, 'FontSize', 16)
ylabel('$\Delta z$ [m]', 'Interpreter', 'LaTeX', 'FontSize', 18)
xlim([0,50])
xlabel('Time step')
ylim(yl)

% Create upper left subplot
subtightplot(2,2,1, gap, marg_h, marg_w)
hold on
plot(kk, zCombinedArray(1:end-1,:),'k-')
hold off
box on
set(gca,'XTick',[], 'FontSize', 14)
ylabel('$\Delta z$ [m]', 'Interpreter', 'LaTeX', 'FontSize', 18)
xlim([0,50])
ylim(yl)
title('$\sigma_{\Delta z} = 3~\rm{km}$', 'Interpreter', 'Latex', 'FontSIze', 20);

% Create lower right subplot
subtightplot(2,2,4,  gap, marg_h, marg_w)
load('./data/run1_arrays.mat');
load('./data/dzArrayPrior.mat');
hold on
priorHandle = plot(0:50, dzArrayPrior, '-', 'Color', 0.7*ones(1,3));
posteriorHandle = plot(0:50, dzArrayPost1, 'k-');
hold off
box on
xlabel('Time step')
legendHandle = legend([priorHandle(1), posteriorHandle(1)], 'Prior', 'Posterior', 'Location', 'NorthWest');
set(gca,'XTick',[0, 10, 20, 30, 40, 50], 'YTick', [],  'FontSize', 14)
set(legendHandle, 'FontSize', 16)
xlim([0,50])
ylim(yl)

% Create upper right subplot
subtightplot(2,2,2, gap, marg_h, marg_w)
hold on
plot(kk, zCombinedArray1(1:end-1,:),'k-')
hold off
box on
set(gca,'XTick',[], 'YTick', [], 'FontSize', 14)
xlim([0,50])
ylim(yl)
title('$\sigma_{\Delta z} = 30~\rm{m}$', 'Interpreter', 'Latex', 'FontSize', 20);

% Reassign variable for next section
grArrayPost = grArray;

%% Smoothed GR observations

dzArrayPost = zwArray - zwArray(1,:);

load('./data/priorResults.mat')
load('./data/Tunalik1grlog.mat')

grArrayPrior = grArray;
dzArrayPrior = [zeros(1,size(dzArray,2)); dzArray];

lightgreen = [0.3, 1.0, 0.3];
lightgreen = [0.0, 0.3, 0.4];
lightred = [1.0, 0.3, 0.3];

gr_smooth = interp1(linspace(0,50,length(gr)), gr, linspace(0,50,512));

grPriorLower = quantile(grArrayPrior',0.10);
grPriorUpper = quantile(grArrayPrior',0.90);

grPostLower = quantile(grArrayPost',0.10);
grPostUpper = quantile(grArrayPost',0.90);

figure()
hold on
ph1 = plot(linspace(0,50,512),gr_smooth,'k-');
fh2 = fill([1:50,fliplr(1:50)], [grPostLower, fliplr(grPostUpper)], lightred);
set(fh2,'EdgeColor','r')
ph2 = plot(grArrayPost(:,4:4:end),'r-');
hold off
ylim([40,100])
box on

legend([ph1, ph2(1)], 'Observed', 'Fitted', 'Location', 'northwest')
xlabel('Time step')
ylabel('GR [API]')
set(gca,'FontSize', 14)

%% Confidence envelopes for thickness estimates

dzPriorLower = quantile(dzArrayPrior',0.10);
dzPriorUpper = quantile(dzArrayPrior',0.90);

dzPostLower = quantile(dzArrayPost',0.10);
dzPostUpper = quantile(dzArrayPost',0.90);

figure()
hold on
fh3 = fill([0:50, fliplr(0:50)], [dzPriorLower, fliplr(dzPriorUpper)], lightgreen);
fh4 = fill([0:50, fliplr(0:50)], [dzPostLower, fliplr(dzPostUpper)], lightred);
set(fh3,'EdgeColor',[0,0.2,0.3])
set(fh4,'EdgeColor','r')
alpha([fh3, fh4], 0.5)
hold off
box on
legend([fh3, fh4], 'A priori', 'A posteriori', 'Location', 'NorthWest')
xlabel('Layer number')
ylabel('Total thickness [m]')
set(gca,'FontSize', 14)
title('Cumulative thickness, 80% confidence')

%% Blockwise GR averages

grBlockMean = [];
for block = 1:10
    grBlockMean(block,:) = mean(grArrayPost((5*block-4):5*block,:),1);
end

figure()
hold on
sth1 = stairs(0:5:50, [grBlockMean; grBlockMean(end,:)], '-', 'Color', 0.8*ones(1,3));
sth2 = stairs(0:5:50, mean([grBlockMean; grBlockMean(end,:)],2), 'k-', 'LineWidth', 2);
sth3 = stairs(0:5:50, mean([grobs;grobs(end,:)],2), 'k--', 'LineWidth', 2);
hold off
legend([sth1(1), sth2, sth3], 'Ensemble members', 'Ensemble mean', 'Observed')
xlabel('Layer')
ylabel('GR [API]')
title('GR log fit, block scale')
box on
set(gca,'FontSize', 14)

figure()
hold on
sth1 = stairs(0:50, [grArrayPost; grArrayPost(end,:)], '-', 'Color', 0.8*ones(1,3));
sth2 = stairs(0:50, mean([grArrayPost; grArrayPost(end,:)],2), 'k-', 'LineWidth', 2);
sth3 = stairs(0:5:50, mean([grobs;grobs(end,:)],2), 'k--', 'LineWidth', 2);
hold off
legend([sth1(1), sth2, sth3], 'Ensemble members', 'Ensemble mean', 'Observed')
xlabel('Layer')
ylabel('GR [API]')
title('GR log fit, layer scale')
box on
set(gca,'FontSize', 14)

load('./data/grArray.mat');

for block = 1:10
    grBlockMean(block,:) = mean(grArray((5*block-4):5*block,:),1);
end

figure()
hold on
sth1 = stairs(0:5:50, [grBlockMean; grBlockMean(end,:)], '-', 'Color', 0.8*ones(1,3));
sth2 = stairs(0:5:50, mean([grBlockMean; grBlockMean(end,:)],2), 'k-', 'LineWidth', 2);
sth3 = stairs(0:5:50, mean([grobs;grobs(end,:)],2), 'k--', 'LineWidth', 2);
hold off
legend([sth1(1), sth2, sth3], 'Prior realizations', 'Prior mean', 'Observed')
xlabel('Layer')
ylabel('GR [API]')
title('GR log prior, block scale')
box on
set(gca,'FontSize', 14)

figure()
hold on
sth1 = stairs(0:50, [grArray; grArray(end,:)], '-', 'Color', 0.8*ones(1,3));
sth2 = stairs(0:50, mean([grArray; grArray(end,:)],2), 'k-', 'LineWidth', 2);
sth3 = stairs(0:5:50, mean([grobs;grobs(end,:)],2), 'k--', 'LineWidth', 2);
hold off
legend([sth1(1), sth2, sth3], 'Prior realizations', 'Prior mean', 'Observed')
xlabel('Layer')
ylabel('GR [API]')
title('GR log prior, layer scale')
box on
set(gca,'FontSize', 14)

%% Vertical plots of synthetic and real GR observations (Fig. 8)

figure()

gap = [0.02, 0.02];
marg_h = [0.02, 0.17];
marg_w = [0.1, 0.05];
ylgr = [25, 118];

subtightplot(1,3,1, gap, marg_h, marg_w)
hold on
sth1 = stairs(0:50, [grArray; grArray(end,:)], '-', 'Color', 0.8*ones(1,3));
sth2 = stairs(0:50, mean([grArray; grArray(end,:)],2), 'k-', 'LineWidth', 2);
sth3 = stairs(0:5:50, mean([grobs;grobs(end,:)],2), 'k--', 'LineWidth', 2);
hold off
legend([sth1(1), sth2, sth3], 'Prior realizations', 'Prior mean', 'Observed', 'Location', 'SouthOutside')
xlabel('Layer')
ylabel('GR [API]')
box on
set(gca,'FontSize', 15, 'YAxisLocation', 'Right')
ylim(ylgr)
view([90, -90])
title({'Prior';'';'';''}, 'Interpreter', 'Latex', 'FontSize', 20)

subtightplot(1,3,2, gap, marg_h, marg_w)
hold on
sth1 = stairs(0:50, [grArrayPost; grArrayPost(end,:)], '-', 'Color', 0.8*ones(1,3));
sth2 = stairs(0:50, mean([grArrayPost; grArrayPost(end,:)],2), 'k-', 'LineWidth', 2);
sth3 = stairs(0:5:50, mean([grobs;grobs(end,:)],2), 'k--', 'LineWidth', 2);
hold off
legend([sth1(1), sth2, sth3], 'Ens. members', 'Ens. mean', 'Observed', 'Location', 'SouthOutside')
ylabel('GR [API]')
box on
set(gca,'FontSize', 15, 'YAxisLocation', 'Right', 'XTick', [])
ylim(ylgr)
view([90, -90])
title({'Posterior';'$\sigma_{\Delta z} = 3~\rm{km}\vspace{1cm}$'}, 'Interpreter', 'Latex', 'FontSize', 20);

subtightplot(1,3,3, gap, marg_h, marg_w)
hold on
sth1 = stairs(0:50, [grArrayPost1; grArrayPost1(end,:)], '-', 'Color', 0.8*ones(1,3));
sth2 = stairs(0:50, mean([grArrayPost1; grArrayPost1(end,:)],2), 'k-', 'LineWidth', 2);
sth3 = stairs(0:5:50, mean([grobs1;grobs1(end,:)],2), 'k--', 'LineWidth', 2);
hold off
legend([sth1(1), sth2, sth3], 'Ens. members', 'Ens. mean', 'Observed', 'Location', 'SouthOutside')
ylabel('GR [API]')
box on
set(gca,'FontSize', 15, 'YAxisLocation', 'Right', 'XTick', [])
ylim(ylgr)
view([90, -90])
title({'Posterior';'$\sigma_{\Delta z} = 30~\rm{m}\vspace{1cm}$'}, 'Interpreter', 'Latex', 'FontSize', 20);

set(gcf, 'Units', 'Centimeters', 'Position', [30,10,25,25])