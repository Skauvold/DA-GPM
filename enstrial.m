% Script for running trials of the Ensemble Smoother (EnS) on the synthetic
% test case described in the paper.


clear; close all; clc

rng(1983)

addpath('../supporting');
addpath('../supporting/jsonlab');

% Number of data assimilation trials to run
numtrials = 10;

% Ensemble properties
ne = 100;  % Ensemble size

% Time axis
tstart = -20000;
tend = 0;
delta_t = 1000;
t = tstart:delta_t:tend;
tlen = length(t);
nsteps = tlen-1;

% Grid properties
nx = 72;
ny = 16;
dims = [nx,ny];

% Additional grid properties
xmin = 0;
xmax = 100*(nx-1);
ymin = 0;
ymax = 100*(ny-1);
xyminmax = [ymin,ymax,xmin,xmax];

% Blind well coordinates
xbw = [10,20,30,40,50,60,70];
ybw = [10,10,10,10,10,10,10];
nbw = length(xbw);

% Initial surface prior parameters
initSurfmean = repmat(linspace(20,-10,nx),ny,1);
bathy_mean = reshape(initSurfmean,prod(dims),1);
sig_z0 = 40;
bathy_cov = sig_z0*sig_z0*ecrm(ny,nx,1,nx/2.5,2);

% Initial sediment proportions prior
isp_prior = 50*[1,1,1,1];

% Sea level parameter prior
sl_prior_mean = linspace(0,30,length(t));
sl_prior_cov = 10*10*ecrm(length(t),1,1,3,2);

% Sediment supply prior parameters
ss_prior_mean = 12*ones(1,length(t));
ss_prior_cov = 5*5*ecrm(length(t),1,1,3,2);
asd_intensity = 0.9999;

% Combined result and reference array
trial_references = cell(numtrials,1);
trial_results = cell(numtrials,1);

% Put prior information into a struct
priorspec = struct;
priorspec.grid.dims = dims;
priorspec.bathy_mean = bathy_mean;
priorspec.bathy_cov = bathy_cov;
priorspec.isp_prior = isp_prior;
priorspec.tstart = tstart;
priorspec.tend = tend;
priorspec.delta_t = delta_t;
priorspec.sl_prior_mean = sl_prior_mean;
priorspec.sl_prior_cov = sl_prior_cov;
priorspec.ss_prior_mean = ss_prior_mean;
priorspec.ss_prior_cov = ss_prior_cov;
priorspec.asd_intensity = asd_intensity;

for trial = 1:numtrials
    % Create reference realization
    reference = mkreffun(priorspec,trial);
    
    % Initial Bathymetry
    Z0ens = mvnrnd(bathy_mean,bathy_cov,ne)'; % ensemble of initial surfaces
    
    % Initial sediment proportions
    p0const = dirchrnd(isp_prior,ne)'; % Global init. sediment proportions
    % Transform to 3-variate real (s).
    s0const = p2eta(p0const);
    
    % Sea level initial ensemble
    f_SL_ens = mvnrnd(sl_prior_mean,sl_prior_cov,ne)';
    
    % Sediment supply initial ensemble
    f_SS_ens = max(0,mvnrnd(ss_prior_mean,ss_prior_cov,ne)');
    area_sed_fun = [asd_intensity*ones(2,ny); zeros(nx-2,ny)];
    
    % Write initial ensemble to files
    for b = 1:ne
        writemap(reshape(Z0ens(:,b),nx,ny),sprintf('input%i_%i',b,trial),p0const(:,b),xyminmax,tstart);
        write_asd_json([f_SS_ens(:,b)';t],area_sed_fun,sprintf('areased%i_%i',b,trial),'../templates/Case5');
        write_flw_json([f_SL_ens(:,b)';t],sprintf('sealevel%i_%i',b,trial),'../templates/Case5');
    end
    
    % Prepare data
    xobs = 50;
    yobs = 6;
    obs_indices = NaN(4,1);
    CY = diag([2,0.5,0.5,0.5]);
    data = struct;
    data.unperturbed = cell(nsteps,1);
    data.perturbed = cell(nsteps,1);
    for i = 1:nsteps
        Y0 = [reference.z{i}((xobs-1)*ny+yobs,end);
            reference.s1{i}((xobs-1)*ny+yobs,end);
            reference.s2{i}((xobs-1)*ny+yobs,end);
            reference.s3{i}((xobs-1)*ny+yobs,end)];
        data.unperturbed{i} = Y0;
        data.perturbed{i} = repmat(Y0,1,ne) + mvnrnd(zeros(1,4),CY,ne)';
    end
    
    % Initialize arrays
    forecast_state = cell(nsteps,1);
    analysis_state = cell(nsteps,1);
    Xf_tmp = cell(1,ne);
    obs_indices_all = [];
    Y = [];
    
    Xa = [f_SL_ens; f_SS_ens; Z0ens; repmat(s0const(1,:),prod(dims),1); repmat(s0const(2,:),prod(dims),1); repmat(s0const(3,:),prod(dims),1)];
    
    tic();
    
    for time_step = 1:nsteps
        nsurfs = time_step;
        Xf = zeros((4*nsurfs+1)*prod(dims)+2*tlen,ne);
        Xf(1:2*tlen,:) = Xa(1:2*tlen,:);
        parfor b = 1:ne
            if time_step >= 2
                f_SL_input = Xa(1:tlen,b);
                f_SS_input = Xa(tlen+1:2*tlen,b);
                modifyOutfile(Xa(2*tlen+1:end,b),sprintf('output%i_%i.out',b,trial),sprintf('input%i_%i.map',b,trial),time_step-1);
                write_asd_json([f_SS_input';t],area_sed_fun,sprintf('areased%i_%i',b,trial),'../templates/Case5');
                write_flw_json([f_SL_input';t],sprintf('sealevel%i_%i',b,trial),'../templates/Case5');
            end
            write_ctl_json([t(time_step),t(time_step+1)],sprintf('enkfcond_%i_%i',b,trial),'../templates/mkref',b,trial);
            % Run simulator [line redacted]
            Xf_tmp{b} = parseOutfile(sprintf('output%i_%i.out',b,trial),time_step);
        end
        for b = 1:ne
            Xf(2*tlen+1:end,b) = Xf_tmp{b};
        end
        forecast_state{time_step} = Xf;
        obs_indices(1) = 2*tlen + nsurfs*prod(dims) + (xobs-1)*ny + yobs; % Zobs
        obs_indices(2) = obs_indices(1) + nsurfs*prod(dims); % S1obs
        obs_indices(3) = obs_indices(2) + nsurfs*prod(dims); % S2obs
        obs_indices(4) = obs_indices(3) + nsurfs*prod(dims); % S3obs
        obs_indices_all = [obs_indices_all;obs_indices];
        
        Y = [Y;data.perturbed{time_step}];
        
        if time_step < nsteps
            Xa = Xf;
        else
            HXf = Xf(obs_indices_all,:);
            HCXfHt = cov(HXf');
            CXfHt = (Xf - repmat(mean(Xf,2),1,ne))*(HXf - repmat(mean(HXf,2),1,ne))'./(ne-1);
            Xa = Xf + CXfHt*((HCXfHt + kron(eye(nsteps),CY))\(Y-HXf));
        end
        analysis_state{time_step} = Xa;
    end
    
    elapsed_time = toc();
    sprintf('Elapsed time: %.2f s\n',elapsed_time)
    
    % Extract relevant part of results
    fresult = struct;
    fresult.bw.xpos = xbw;
    fresult.bw.ypos = ybw;
    fresult.z = cell(tlen-1,nbw);
    fresult.s1 = cell(tlen-1,nbw);
    fresult.s2 = cell(tlen-1,nbw);
    fresult.s3 = cell(tlen-1,nbw);
    fresult.thetaSL = cell(1,tlen);
    fresult.thetaSS = cell(1,tlen);
    
    smallref = struct;
    smallref.z = cell(tlen-1,nbw);
    smallref.s1 = cell(tlen-1,nbw);
    smallref.s2 = cell(tlen-1,nbw);
    smallref.s3 = cell(tlen-1,nbw);
    smallref.thetaSL = reference.seafun;
    smallref.thetaSS = reference.sedfun;
    
    for k = 1:nsteps
        for ibw = 1:nbw
            bwinds_z = 2*tlen + (xbw(ibw)-1)*ny + ybw(ibw) + prod(dims)*[0:k];
            bwinds_s1 = 2*tlen + (k+1)*prod(dims) + (xbw(ibw)-1)*ny + ybw(ibw) + prod(dims)*[0:(k-1)];
            bwinds_s2 = 2*tlen + (2*k+1)*prod(dims) + (xbw(ibw)-1)*ny + ybw(ibw) + prod(dims)*[0:(k-1)];
            bwinds_s3 = 2*tlen + (3*k+1)*prod(dims) + (xbw(ibw)-1)*ny + ybw(ibw) + prod(dims)*[0:(k-1)];
            fresult.z{k,ibw} = analysis_state{k}(bwinds_z,:);
            fresult.s1{k,ibw} = analysis_state{k}(bwinds_s1,:);
            fresult.s2{k,ibw} = analysis_state{k}(bwinds_s2,:);
            fresult.s3{k,ibw} = analysis_state{k}(bwinds_s3,:);
            
            bwinds_ref = (xbw(ibw)-1)*ny + ybw(ibw);
            smallref.z{k,ibw} = reference.z{k}(bwinds_ref,:)';
            smallref.s1{k,ibw} = reference.s1{k}(bwinds_ref,:)';
            smallref.s2{k,ibw} = reference.s2{k}(bwinds_ref,:)';
            smallref.s3{k,ibw} = reference.s3{k}(bwinds_ref,:)';
        end
        fresult.thetaSL{k} = analysis_state{k}(1:tlen,:);
        fresult.thetaSS{k} = analysis_state{k}(tlen+1:2*tlen,:);
    end
    
    trial_references{trial} = smallref;
    trial_results{trial} = fresult;
    
    % Delete almost everything
    delete *.map
    delete *.out
    delete *.flw.json
    delete *.asd.json
    delete *.ctl.json
    
    % Copy dfu and sed files back into working dir.
    copyfile('../templates/mkref.dfu.json','.');
    copyfile('../templates/mkref.sed.json','.');
    
    % Store results of individual trial
    save_partial_results(sprintf('results_trial%i',trial),smallref,fresult);
end

% Store final results
save('trial_output','trial_references','trial_results');

%cleanup
delete *.map
delete *.out
delete *.flw.json
delete *.asd.json
delete *.ctl.json