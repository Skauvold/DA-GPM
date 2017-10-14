% Use the EnKF to condition North Slope, Alaska simulation to gamma ray
% data from the Tunalik 1 well.

clear all; close all; clc
rng(56784928)   % Set seed for random number generator

% Add required directories to path
addpath('../supporting');
addpath('../supporting/jsonlab');
addpath('../supporting/Inpaint_nans');

% Load observations and initial condition data from files
load('../data/ANSdata_Tunalik1_20170116.mat')
load('../data/tunalik1InputData.mat')

% Grid dimensions
dims = reference.grid.dims;
nx = dims(1);
ny = dims(2);
xmin = 0;
xmax = 15000*(nx-1);
ymin = 0;
ymax = 15000*(ny-1);
xyminmax = [ymin,ymax,xmin,xmax];

% Time axis
tstart = reference.time.start;
tend = reference.time.end;
delta_t = reference.time.steplen;
t = tstart:delta_t:tend;
nsteps = reference.time.nsteps;
tlen = length(t);

% Specify when to update the ensemble
conditioning_time_steps = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50];

% Ensemble size
ne = 100;

% Observation properties
xobs = reference.well.xind - 3;
yobs = reference.well.yind - 2;
obsSd = [30, 3, 0.005]; % (z, gr_average, gr_slope)
grValues = [25, 40, 55, 115]; % Representative GR values

% Initial Bathymetry
initSurfmean = reference.initialSurface; % Mean surface
zoffset = mvnrnd([0,0,0],diag(120*120*ones(1,3)),ne)'; % Perturbations
Z0ens = NaN(prod(dims),ne);

% Create perturbed surfaces
for b = 1:ne
    [msgr1, msgr2] = meshgrid(linspace(zoffset(1,b),zoffset(2,b),ny), ...
        linspace(zoffset(1,b),zoffset(3,b),nx));
    z0b = initSurfmean + msgr1 + msgr2;
    Z0ens(:,b) = reshape(z0b',prod(dims),1);
end

% Initial sediment proportions
p0const = dirchrnd(0.5*[10,5,1,1],ne)'; % Global init. sediment proportions
% Transform to 3-variate real (s).
s0const = p2eta(p0const);

% Sea level initial ensemble
f_SL_ens = NaN(length(t),ne);
for b = 1:ne
    midb = randi(round(length(t).*[0.3,0.5]));
    SLb_start = normrnd(110,1400);
    SLb_end = normrnd(160,1200);
    SLb_mid = normrnd(0.5*(SLb_start+SLb_end),800);
    fSLb = [linspace(SLb_start,SLb_mid,midb+1)';
                     linspace(SLb_mid,SLb_end,length(t)-midb)']-350;
    fSLb(midb+1) = [];
    f_SL_ens(:,b) = fSLb;
end

% Specifty grid properties
xcorners = [-1200000, 435000, 435000, -1200000];
ycorners = [1400000, 1400000, 2690000, 2690000];
corners = [xcorners;ycorners];

% Define source map
SourceLocations = [
82, 25;
78, 25;
74, 25;
70, 25;
65, 25;
58, 26;
53, 26;
50, 27;
44, 28;
40, 30;
36, 30;
32, 32;
27, 36;
22, 41;
20, 46;
18, 53;
17, 60;
17, 64;
17, 68;
16, 72;
16, 78;
16, 85;
16, 89;
16, 94;
16, 100;
16, 105;
];
mySourceMap = makeSourceMap([nx,ny],SourceLocations,3,0.9999);
imagesc(mySourceMap)
zeroSourceMap = zeros(nx,ny);

% Sediment supply parameter
f_SS_ens = NaN(length(t),ne);
for b = 1:ne
    f_SS_start = min(abs(normrnd(0.4, 0.2)),0.67);
    f_SS_end = min(abs(normrnd(0.4, 0.2)),0.67);
    f_SS_ens(:,b) = linspace(f_SS_start,f_SS_end,length(t))';
end

% Write initial ensemble to files
for b = 1:ne
    writemapAndSrc(reshape(Z0ens(:,b),nx,ny),sprintf('input%i',b), ...
        p0const(:,b),xyminmax,tstart,zeroSourceMap,corners);
    write_flw_json([f_SL_ens(:,b)';t],sprintf('sealevel%i',b),'../templates/ANSCase');
    write_asd_json([f_SS_ens(:,b)';t],mySourceMap,sprintf('areased%i',b),'../templates/NSAKTemplate');
end

% Initialize array
saveIndexVector = conditioning_time_steps;
saveIndexCounter = 0;
synthetic_apriori = cell(length(saveIndexVector),1);
synthetic_aposteriori = cell(length(saveIndexVector),1);
observed_data = cell(length(saveIndexVector),1);

Xa = [f_SL_ens; f_SS_ens; Z0ens; repmat(s0const(1,:),prod(dims),1); repmat(s0const(2,:),prod(dims),1); repmat(s0const(3,:),prod(dims),1)];

tic();

for time_step = 1:nsteps
    nsurfs = time_step;
    Xf = zeros((4*nsurfs+1)*prod(dims)+2*tlen,ne);
    Xf(1:2*tlen,:) = Xa(1:2*tlen,:); % Copy parameters (SL, SS) to Xf
    parfor b = 1:ne
        if time_step >= 2
            f_SL_input = Xa(1:tlen,b);
            f_SS_input = Xa((tlen+1):2*tlen,b);
            modifyOutfile4(Xa(2*tlen+1:end,b),sprintf('output%i.out',b),sprintf('input%i.map',b),time_step-1);
            write_flw_json([f_SL_input';t],sprintf('sealevel%i',b),'../templates/ANSCase');
            write_asd_json([f_SS_input';t],mySourceMap,sprintf('areased%i',b),'../templates/NSAKTemplate');
        end
        write_ctl_json([t(time_step),t(time_step+1)],sprintf('ansenkf_%i',b),'../templates/ANSCase',b);
        % Run simulator [line redacted]
        Xft{b} = parseOutfile4(sprintf('output%i.out',b),time_step);
    end
    for b = 1:ne
        Xf(2*tlen+1:end,b) = Xft{b};
    end
    fprintf('Forecast step %i of %i complete.\n\n',time_step,nsteps)
    if any(time_step == conditioning_time_steps)
        cts_index = find(conditioning_time_steps == time_step);
        if cts_index == 1
            condInterval = [1:time_step];
        else
            condInterval = [(conditioning_time_steps(cts_index-1)+1):time_step];
        end        
        HXf = simplyGetHx(Xf, condInterval, tlen, dims, grValues, [xobs, yobs]);
        [Y, CY] = simplyGetY(tunalik1InputData, condInterval, obsSd, ne);
        HCXfHt = cov(HXf');
        CXfHt = (Xf - repmat(mean(Xf,2),1,ne))*(HXf - repmat(mean(HXf,2),1,ne))'./(ne-1);
        Xa = Xf + CXfHt*((HCXfHt + CY)\(Y - HXf));
        Xa(tlen+1:2*tlen,:) = min(max(Xa(tlen+1:2*tlen,:),zeros(tlen,ne)),0.67*ones(tlen,ne));
        fprintf('Layers %i-%i of %i updated.\n\n\n',condInterval(1),condInterval(end),nsteps)
    else
        Xa = Xf;
    end
    if any(saveIndexVector == time_step) % Store results at given times
        saveIndexCounter = saveIndexCounter + 1;
        synthetic_apriori{saveIndexCounter} = HXf;
        HXa = simplyGetHx(Xa, condInterval, tlen, dims, grValues, [xobs, yobs]);
        synthetic_aposteriori{saveIndexCounter} = HXa;
        observed_data{saveIndexCounter} = Y;
    end
end    

elapsed_time = toc()

% Store final results
save('ANSCaseEnKF_results','observed_data','synthetic_apriori','synthetic_aposteriori','Xa');

% Cleanup working directory
delete *.map
delete *.out
delete *.flw.json
delete *.asd.json
delete *.ctl.json
