function reference = mkreffun(spec,idnum);
% A function for generating a reference realization by running GPM

% Grid properties
dims = spec.grid.dims;
nx = dims(1);
ny = dims(2);
xx = 1:nx;
yy = 1:ny;

xmin = 0;
xmax = 100*(nx-1);
ymin = 0;
ymax = 100*(ny-1);

xyminmax = [ymin,ymax,xmin,xmax];

% Initial Bathymetry
z0vec = mvnrnd(spec.bathy_mean,spec.bathy_cov);

% Initial sediment proportions
p0vec = dirchrnd(spec.isp_prior,1)';
% Transform to 3-variate real (s).
s0vec = p2eta(p0vec);

% Time axis and number of time steps
tstart = spec.tstart;
tend = spec.tend;
delta_t = spec.delta_t;
t = tstart:delta_t:tend;
nsteps = length(t)-1;

% Sea level
f_SL = mvnrnd(spec.sl_prior_mean,spec.sl_prior_cov);
write_flw_json([f_SL;t],sprintf('sealevel%i',idnum),'../templates/Case5');

% Sediment supply
f_SS = max(0,mvnrnd(spec.ss_prior_mean,spec.ss_prior_cov));
area_sed_fun = [spec.asd_intensity*ones(2,ny); zeros(nx-2,ny)];   
write_asd_json([f_SS;t],area_sed_fun,sprintf('areased%i',idnum),'../templates/Case5');

% Simulation

% Forecast state cell array
Xf = cell(nsteps,1);

% Prepare mapfile
writemap(reshape(z0vec,nx,ny),sprintf('input%i',idnum),p0vec,xyminmax,tstart);

for k = 1:nsteps
    write_ctl_json([t(k),t(k+1)],sprintf('mkref%i',idnum),'../templates/mkref',idnum);
    % Run simulator [line redacted]
    Xf{k} = parseOutfile(sprintf('output%i.out',idnum),k);
    movefile(sprintf('output%i.out',idnum),sprintf('input%i.map',idnum));
end

% Organize reference data as a struct
reference = struct;
reference.p1 = cell(1,nsteps);
reference.p2 = cell(1,nsteps);
reference.p3 = cell(1,nsteps);
reference.s1 = cell(1,nsteps);
reference.s2 = cell(1,nsteps);
reference.s3 = cell(1,nsteps);
reference.z = cell(1,nsteps);
reference.seafun = f_SL;
reference.sedfun = f_SS;
reference.time = struct;
reference.time.start = tstart;
reference.time.end = tend;
reference.time.steplen = delta_t;
reference.time.nsteps = nsteps;
reference.grid = struct;
reference.grid.dims = dims;
reference.grid.xlims = [xmin,xmax];
reference.grid.ylims = [ymin,ymax];

% Store reference realisation in struct
for i = 1:nsteps
    s1_i = reshape(Xf{i}((i+1)*prod(dims)+1:(2*i+1)*prod(dims)),prod(dims),i);
    s2_i = reshape(Xf{i}((2*i+1)*prod(dims)+1:(3*i+1)*prod(dims)),prod(dims),i);
    s3_i = reshape(Xf{i}((3*i+1)*prod(dims)+1:(4*i+1)*prod(dims)),prod(dims),i);
    reference.s1{i} = s1_i;
    reference.s2{i} = s2_i;
    reference.s3{i} = s3_i;
    p_i = eta2p([Xf{i}((i+1)*prod(dims)+1:(2*i+1)*prod(dims))'; ...
                 Xf{i}((2*i+1)*prod(dims)+1:(3*i+1)*prod(dims))';
                 Xf{i}((3*i+1)*prod(dims)+1:(4*i+1)*prod(dims))']);
    reference.p1{i} = reshape(p_i(1,:),prod(dims),i)';
    reference.p2{i} = reshape(p_i(2,:),prod(dims),i)';
    reference.p3{i} = reshape(p_i(3,:),prod(dims),i)';
    reference.z{i} = reshape(Xf{i}(1:(i+1)*prod(dims)),prod(dims),i+1);
end

end
