% Compute summary statistics for EnKF or EnS trial results

clear; close all; clc;

Ntr = 177; % Number of trials
Nw = 7; % Number of "wells"
nt = 20; % Number of layers
B = 100; % Number of ensemble members
alpha_ci = 0.2; % Confidence interval significance level

% setup arrays
sqeZ = NaN(nt,Nw,Ntr);
sqeS = NaN(nt,Nw,Ntr);

crpsZ = NaN(nt,Nw,Ntr);
crpsS1 = NaN(nt,Nw,Ntr);
crpsS2 = NaN(nt,Nw,Ntr);
crpsS3 = NaN(nt,Nw,Ntr);

cvgZ = NaN(nt,Nw,Ntr);
cvgS1 = NaN(nt,Nw,Ntr);
cvgS2 = NaN(nt,Nw,Ntr);
cvgS3 = NaN(nt,Nw,Ntr);

sqeSL = NaN(nt+1,Ntr);
sqeSS = NaN(nt+1,Ntr);

crpsSL = NaN(nt+1,Ntr);
crpsSS = NaN(nt+1,Ntr);

cvgSL = NaN(nt+1,Ntr);
cvgSS = NaN(nt+1,Ntr);

tmark = tic;

% Step through trial result files and extract required data
for tr = 1:Ntr
    load(sprintf('../Sixth_run/results_trial%i.mat',tr));
    for w = 1:Nw
        for k = 1:nt
            % Squared error
            sqeZ(k,w,tr) = (mean(res.z{nt,w}(k,:)) - ref.z{nt,w}(k))^2;
            sqeS1 = (mean(res.s1{nt,w}(k,:)) - ref.s1{nt,w}(k))^2;
            sqeS2 = (mean(res.s2{nt,w}(k,:)) - ref.s2{nt,w}(k))^2;
            sqeS3 = (mean(res.s3{nt,w}(k,:)) - ref.s3{nt,w}(k))^2;
            sqeS(k,w,tr) = sqeS1 + sqeS2 + sqeS3;
            
            % CRPS
            crpsZ(k,w,tr) = crps(res.z{nt,w}(k,:),ref.z{nt,w}(k));
            % How to compute CRPS for s = (s1,s2,s3)?
            % CI or confidence ellipsoid for s?
            crpsS1(k,w,tr) = crps(res.s1{nt,w}(k,:),ref.s1{nt,w}(k));
            crpsS2(k,w,tr) = crps(res.s2{nt,w}(k,:),ref.s2{nt,w}(k));
            crpsS3(k,w,tr) = crps(res.s3{nt,w}(k,:),ref.s3{nt,w}(k));
            
            % CI coverage of Z
            Lz = quantile(res.z{nt,w}(k,:),alpha_ci/2);
            Uz = quantile(res.z{nt,w}(k,:),1-alpha_ci/2);
            Tz = ref.z{nt,w}(k);
            if Lz <= Tz && Tz <= Uz
                cvgZ(k,w,tr) = 1;
            else
                cvgZ(k,w,tr) = 0;
            end
            
            % CI coverage of S1
            LS1 = quantile(res.s1{nt,w}(k,:),alpha_ci/2);
            US1 = quantile(res.s1{nt,w}(k,:),1-alpha_ci/2);
            TS1 = ref.s1{nt,w}(k);
            if LS1 <= TS1 && TS1 <= US1
                cvgS1(k,w,tr) = 1;
            else
                cvgS1(k,w,tr) = 0;
            end
            
            % CI coverage of S2
            LS2 = quantile(res.s2{nt,w}(k,:),alpha_ci/2);
            US2 = quantile(res.s2{nt,w}(k,:),1-alpha_ci/2);
            TS2 = ref.s2{nt,w}(k);
            if LS2 <= TS2 && TS2 <= US2
                cvgS2(k,w,tr) = 1;
            else
                cvgS2(k,w,tr) = 0;
            end
            
            % CI coverage of S3
            LS3 = quantile(res.s3{nt,w}(k,:),alpha_ci/2);
            US3 = quantile(res.s3{nt,w}(k,:),1-alpha_ci/2);
            TS3 = ref.s3{nt,w}(k);
            if LS3 <= TS3 && TS3 <= US3
                cvgS3(k,w,tr) = 1;
            else
                cvgS3(k,w,tr) = 0;
            end
        end
    end
    
    % Compute summary statistics for single trial
    
    % squared error for parameters
    sqeSL(:,tr) = (mean(res.thetaSL{nt},2) - ref.thetaSL').^2;
    sqeSS(:,tr) = (mean(res.thetaSS{nt},2) - ref.thetaSS').^2;
    
    % CRPS for parameters
    for k = 1:nt+1
        crpsSL(k,tr) = crps(res.thetaSL{nt}(k,:),ref.thetaSL(k));
        crpsSS(k,tr) = crps(res.thetaSS{nt}(k,:),ref.thetaSS(k));
    end
    
    % CI coverage for parameters
    LSL = quantile(res.thetaSL{nt},alpha_ci/2,2);
    USL = quantile(res.thetaSL{nt},1-alpha_ci/2,2);
    TSL = ref.thetaSL';
    cvgSL(:,tr) = LSL <= TSL & TSL <= USL;
    
    LSS = quantile(res.thetaSS{nt},alpha_ci/2,2);
    USS = quantile(res.thetaSS{nt},1-alpha_ci/2,2);
    TSS = ref.thetaSS';
    cvgSS(:,tr) = LSS <= TSS & TSS <= USS;
    
    LSLarray(:,tr) = LSL;
    USLarray(:,tr) = USL;
    TSLarray(:,tr) = TSL;
    
    LSSarray(:,tr) = LSS;
    USSarray(:,tr) = USS;
    TSSarray(:,tr) = TSS;
end

elapsed = toc(tmark);
sprintf('Elapsed time: %.2f s\n',elapsed)

%% Organise data

% Discard spurious realizations
keeps0 = squeeze(log10(max(range(sqeZ)))<10);
keeps = find(keeps0==1);
keeps = keeps(1:100);

sqeZ = sqeZ(:,:,keeps);
sqeS = sqeS(:,:,keeps);

crpsZ = crpsZ(:,:,keeps);
crpsS1 = crpsS1(:,:,keeps);
crpsS2 = crpsS2(:,:,keeps);
crpsS3 = crpsS3(:,:,keeps);

cvgZ = cvgZ(:,:,keeps);
cvgS1 = cvgS1(:,:,keeps);
cvgS2 = cvgS2(:,:,keeps);
cvgS3 = cvgS3(:,:,keeps);

sqeSL = sqeSL(:,keeps);
sqeSS = sqeSS(:,keeps);

crpsSL = crpsSL(:,keeps);
crpsSS = crpsSS(:,keeps);

cvgSL = cvgSL(:,keeps);
cvgSS = cvgSS(:,keeps);

% Compute aggregate summary statistics for all trials

% Compute MSE, average CRPS and coverage probability

mseZ = mean(sqeZ,3);
mseS = mean(sqeS,3);

mcrpsZ = mean(crpsZ,3);
mcrpsS1 = mean(crpsS1,3);
mcrpsS2 = mean(crpsS2,3);
mcrpsS3 = mean(crpsS3,3);

cvprZ = mean(cvgZ,3);
cvprS1 = mean(cvgS1,3);
cvprS2 = mean(cvgS2,3);
cvprS3 = mean(cvgS3,3);

mseSL = mean(sqeSL,2);
mseSS = mean(sqeSS,2);

mcrpsSL = mean(crpsSL,2);
mcrpsSS = mean(crpsSS,2);

cvprSL = mean(cvgSL,2);
cvprSS = mean(cvgSS,2);

USLarray = USLarray(:, keeps);
TSLarray = TSLarray(:, keeps);
LSLarray = LSLarray(:, keeps);

USSarray = USSarray(:, keeps);
TSSarray = TSSarray(:, keeps);
LSSarray = LSSarray(:, keeps);

%% Make table of well-wise MSE, CRPS and Cov. Prob. for Z, S, thetaSL and thetaSS
clc

% Z
fprintf('& %.2f\t',mean(mseZ,1))
fprintf('\\\\ \n')
fprintf('& %.2f\t',mean(mcrpsZ,1))
fprintf('\\\\ \n')
fprintf('& %.2f\t', mean(cvprZ,1))
fprintf('\\\\ \n\n')

% S
mcrpsS = (mcrpsS1 + mcrpsS2 + mcrpsS3)/3;
cvprS = (cvprS1 + cvprS2 + cvprS3)/3;
fprintf('& %.2f\t',mean(mseS,1))
fprintf('\\\\ \n')
fprintf('& %.2f\t',mean(mcrpsS,1))
fprintf('\\\\ \n')
fprintf('& %.2f\t', mean(cvprS,1))
fprintf('\\\\ \n\n')

%% Make tables of MSE, CRPS and Cov. Prob. for Z, S, thetaSL and thetaSS
clc

% Z
fprintf('%.2f\n',mean(mean(mseZ)))
fprintf('%.2f\n',mean(mean(mcrpsZ)))
fprintf('%.2f\n\n', mean(mean(cvprZ)))

mcrpsS = (mcrpsS1 + mcrpsS2 + mcrpsS3)/3;
cvprS = (cvprS1 + cvprS2 + cvprS3)/3;

% S
fprintf('%.2f\n',mean(mean(mseS)))
fprintf('%.2f\n',mean(mean(mcrpsS)))
fprintf('%.2f\n\n', mean(mean(cvprS)))

% theta_SL
fprintf('%.2f\n',mean(mean(mseSL)))
fprintf('%.2f\n',mean(mean(mcrpsSL)))
fprintf('%.2f\n\n', mean(mean(cvprSL)))

% theta_SS
fprintf('%.2f\n',mean(mean(mseSS)))
fprintf('%.2f\n',mean(mean(mcrpsSS)))
fprintf('%.2f\n\n', mean(mean(cvprSS)))

%% More tables
clc

fprintf('MSE for Z\n')
fprintf('\\hline\n')
fprintf('\t& %i',1:7)
fprintf('\\\\ \n')
fprintf('\\hline\n')
for k = nt:-1:1
    fprintf('%i',k)
    fprintf('\t& $%.2f$',mseZ(k,:))
    fprintf(' \\\\ \n')
end
fprintf('\\hline\n')

fprintf('\nMSE for S\n')
fprintf('\\hline\n')
fprintf('\t& %i',1:7)
fprintf('\\\\ \n')
fprintf('\\hline\n')
for k = nt:-1:1
    fprintf('%i',k)
    fprintf('\t& $%.2f$',mseS(k,:))
    fprintf(' \\\\ \n')
end
fprintf('\\hline\n')

fprintf('\nAvg. CRPS for Z\n')
fprintf('\\hline\n')
fprintf('\t& %i',1:7)
fprintf('\\\\ \n')
fprintf('\\hline\n')
for k = nt:-1:1
    fprintf('%i',k)
    fprintf('\t& $%.2f$',mcrpsZ(k,:))
    fprintf(' \\\\ \n')
end
fprintf('\\hline\n')

fprintf('\nCVPR for Z\n')
fprintf('\\hline\n')
fprintf('\t& %i',1:7)
fprintf('\\\\ \n')
fprintf('\\hline\n')
for k = nt:-1:1
    fprintf('%i',k)
    fprintf('\t& $%.2f$',cvprZ(k,:))
    fprintf(' \\\\ \n')
end
fprintf('\\hline\n')

fprintf('\nMSE/CRPS/CVPR for SL/SS\n')
fprintf('\\hline \n')
fprintf('&\tMSE\t& CRPS\t& Cov. pr. &\tMSE\t& CRPS\t& Cov. pr. \\\\\n')
fprintf('\\hline \n')
for k = 1:nt
    fprintf('%i',k)
    fprintf('\t& $%.2f$',[mseSL(k),mcrpsSL(k),cvprSL(k),mseSS(k),mcrpsSS(k),cvprSS(k)])
    fprintf(' \\\\ \n')
end
fprintf('\\hline \n')

%% Make images displaying CRPS of S1, S2, S3

fh = figure();
ax(1) = subplot(1,3,1);
imagesc(mcrpsS1,[0,10])
set(gca,'ydir','normal','FontSize',14)
ylabel('Layer number')
title('$S_1$','Interpreter','Latex')

ax(2) = subplot(1,3,2);
imagesc(mcrpsS2,[0,10])
set(gca,'ydir','normal','ytick',[],'Yticklabel',[],'FontSize',14)
xlabel('Well number')
title('$S_2$','Interpreter','Latex')

ax(3) = subplot(1,3,3);
imagesc(mcrpsS3,[0,10])
set(gca,'ydir','normal','Ytick',[],'Yticklabel',[],'FontSize',14)
title('$S_3$','Interpreter','Latex')

ch=colorbar;
ylabel(ch,'Average CRPS','FontSize',16)
set(ch, 'Position', [.82 .11 .05 .8150])
for i=1:3
      pos=get(ax(i), 'Position');
      set(ax(i), 'Position', [0.85*pos(1) pos(2) 0.95*pos(3) pos(4)]);
end

colormap(gray);
set(ch,'FontSize',14);

%% Make images displaying coverage probability of S1, S2, S3

fh = figure();
ax(1) = subplot(1,3,1);
imagesc(cvprS1,[0,1])
set(gca,'ydir','normal','FontSize',14)
ylabel('Layer number')
title('$S_1$','Interpreter','Latex')

ax(2) = subplot(1,3,2);
imagesc(cvprS2,[0,1])
set(gca,'ydir','normal','ytick',[],'Yticklabel',[],'FontSize',14)
xlabel('Well number')
title('$S_2$','Interpreter','Latex')

ax(3) = subplot(1,3,3);
imagesc(cvprS3,[0,1])
set(gca,'ydir','normal','Ytick',[],'Yticklabel',[],'FontSize',14)
title('$S_3$','Interpreter','Latex')

ch=colorbar;
ylabel(ch,'Coverage probability','FontSize',16)
set(ch, 'Position', [.82 .11 .05 .8150])
for i=1:3
      pos=get(ax(i), 'Position');
      set(ax(i), 'Position', [0.85*pos(1) pos(2) 0.95*pos(3) pos(4)]);
end

colormap(gray);
set(ch,'FontSize',14);

%% Make figures of theta_SL and theta_SS

tt = linspace(-20,0,nt+1);

figure();

load('LUTSSLL_EnKF.mat')
subplot(2,2,1)
hold on
filh = fill([tt,fliplr(tt)],[LSS;flipud(USS)],0.8*[1,1,1]);
set(filh,'edgecolor','none')
plot(tt,TSS,'k-','LineWidth',1);
title('EnKF');
hold off
box on
ylabel('$\theta_\mathrm{SS}$','Interpreter','Latex','FontSize',16);
set(gca,'XtickLabel',[],'FontSize',14);

subplot(2,2,3)
hold on
filh = fill([tt,fliplr(tt)],[LSL;flipud(USL)],0.8*[1,1,1]);
set(filh,'edgecolor','none')
plot(tt,TSL,'k-','LineWidth',1);
hold off
box on
xlabel('Time (ky)','Interpreter','Latex','FontSize',14);
ylabel('$\theta_\mathrm{SL}$','Interpreter','Latex','FontSize',16);
set(gca,'FontSize',14)

load('LUTSSLL_EnS.mat')
subplot(2,2,2)
hold on
filh = fill([tt,fliplr(tt)],[LSS;flipud(USS)],0.8*[1,1,1]);
set(filh,'edgecolor','none')
plot(tt,TSS,'k-','LineWidth',1);
title('EnS');
hold off
box on
ylabel('$\theta_\mathrm{SS}$','Interpreter','Latex','FontSize',16);
set(gca,'XtickLabel',[],'FontSize',14);

subplot(2,2,4)
hold on
filh = fill([tt,fliplr(tt)],[LSL;flipud(USL)],0.8*[1,1,1]);
set(filh,'edgecolor','none')
plot(tt,TSL,'k-','LineWidth',1);
hold off
box on
xlabel('Time (ky)','Interpreter','Latex','FontSize',14);
ylabel('$\theta_\mathrm{SL}$','Interpreter','Latex','FontSize',16);
set(gca,'FontSize',14);

%% Make plots of MSE/CRPS/CVPR for SL and SS

load('SLSS_stats.mat')
load('SLSS_stats_EnS.mat')

figure()
subplot(3,2,1)
hold on
plot(tt,mseSS_EnKF,'k-')
plot(tt,mseSS,'k--')
hold off
ylabel('$\mathrm{MSE}(\theta_\mathrm{SS})$','Interpreter','Latex','FontSize',14)
set(gca,'XTickLabel',[],'FontSize',14)
title('Sediment Supply')

subplot(3,2,3)
hold on
plot(tt,mcrpsSS_EnKF,'k-')
plot(tt,mcrpsSS,'k--')
hold off
ylabel('$\mathrm{CRPS}(\theta_\mathrm{SS})$','Interpreter','Latex','FontSize',14)
set(gca,'XTickLabel',[],'FontSize',14)
ylim([0,20])

subplot(3,2,5)
hold on
plot(tt,cvprSS_EnKF,'k-')
plot(tt,cvprSS,'k--')
hold off
ylabel('$\mathrm{C.pr.}(\theta_\mathrm{SS})$','Interpreter','Latex','FontSize',14)
xlabel('Time (ky)','Interpreter','Latex','FontSize',14)
set(gca,'FontSize',14)
ylim([0,1])

subplot(3,2,2)
hold on
plot(tt,mseSL_EnKF,'k-')
plot(tt,mseSL,'k--')
hold off
ylabel('$\mathrm{MSE}(\theta_\mathrm{SL})$','Interpreter','Latex','FontSize',14)
set(gca,'XTickLabel',[],'FontSize',14)
title('Sea Level')

subplot(3,2,4)
hold on
plot(tt,mcrpsSL_EnKF,'k-')
plot(tt,mcrpsSL,'k--')
hold off
ylabel('$\mathrm{CRPS}(\theta_\mathrm{SL})$','Interpreter','Latex','FontSize',14)
set(gca,'XTickLabel',[],'FontSize',14)
ylim([0,20])

subplot(3,2,6)
hold on
plot(tt,cvprSL_EnKF,'k-')
plot(tt,cvprSL,'k--')
hold off
ylabel('$\mathrm{C.pr.}(\theta_\mathrm{SL})$','Interpreter','Latex','FontSize',14)
xlabel('Time (ky)','Interpreter','Latex','FontSize',14)
set(gca,'FontSize',14)
ylim([0,1])