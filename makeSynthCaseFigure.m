% A script for creating Fig. 2

clear all; close all; clc

% Load reference realisation
load('./data/myref.mat');

% Set grid properties
k = 0;
nsteps = reference.time.nsteps;
dims = reference.grid.dims;
nx = dims(1);
ny = dims(2);

% Define colours
Color1 = [255, 0, 0]./255; % red
Color2 = [255, 246, 0]./255; % yellow
Color3 = [60, 206, 30]./255; % green
ColorMatrix = [Color1; Color2; Color3];

% Set time step at which to plot layer package
step = 20;

% Create figure window
figure()
hold on

% Step through layers and plot external cell faces
for i = 2:20
    currentSurf = reshape(reference.z{step}(:,i),dims(2),dims(1));
    nextlowerSurf = reshape(reference.z{step}(:,i-1),dims(2),dims(1));
    currentGrain1 = reshape(reference.p1{step}(i-1,:),dims(2),dims(1));
    currentGrain2 = reshape(reference.p2{step}(i-1,:),dims(2),dims(1));
    currentGrain3 = reshape(reference.p3{step}(i-1,:),dims(2),dims(1));
    for j = 2:dims(1)
        xj = [j-1, j, j, j-1];
        yj = [1, 1, 1, 1];
        zj = [nextlowerSurf(1,j-1:j),fliplr(currentSurf(1,j-1:j))];
        PjSand = currentGrain1(1,j-1) + currentGrain2(1,j-1);
        PjSilt = currentGrain3(1,j-1);
        PjClay = 1 - PjSand - PjSilt;
        Pj = [PjSand, PjSilt, PjClay];
        cj = Pj*ColorMatrix;
        fill3(xj,yj,zj,cj);
    end
    for k = 2:dims(2)
        xk = [nx, nx, nx, nx];
        yk = [k-1, k, k, k-1];
        zk = [nextlowerSurf(k-1:k,nx)',fliplr(currentSurf(k-1:k,nx)')];
        PkSand = currentGrain1(k-1,nx) + currentGrain2(k-1,nx);
        PkSilt = currentGrain3(k-1,nx);
        PkClay = 1 - PkSand - PkSilt;
        Pk = [PkSand, PkSilt, PkClay];
        ck = Pk*ColorMatrix;
        fill3(xk,yk,zk,ck);
    end
    if i==nsteps
        for j = 2:dims(1)
            for k = 2:dims(2)
                xjk = [j-1, j, j, j-1];
                yjk = [k-1, k-1, k, k];
                zjk = [currentSurf(k-1,j-1), currentSurf(k-1,j), currentSurf(k,j), currentSurf(k,j-1)];
                PjkSand = currentGrain1(k,j) + currentGrain2(k,j);
                PjkSilt = currentGrain3(k,j);
                PjkClay = 1 - PjkSand - PjkSilt;
                Pjk = [PjkSand, PjkSilt, PjkClay];
                cjk = Pjk*ColorMatrix;
                fill3(xjk, yjk, zjk, cjk);
            end
        end
    end
end
view(20,3)
xlim([3,dims(1)])
ylim([0,dims(2)])
zlim([-20,60])

% Set aspect ratio of axes
daspect([1,1,2])
view(40, 15);

% Add ticks to axes
set(gca, 'XTick', 10:10:70, 'XTickLabel', 1:7)
set(gca, 'YTick', [0, 15], 'YTickLabel', [0, 1.5])

% Add labels and set font size
xlabel('x [km]')
ylabel('y [km]')
zlabel('Elevation [m]')
set(gca, 'FontSize', 16)

%% Make legend

% Create a new set of axes in the figure
axes('Position',[.65 .65 .2 .2])
set(gca, 'visible', 'off')

% Define coordinates for legend vertices
myXlims = get(gca,'Xlim');
xtr = myXlims(1)*ones(1,3) + [0.6, 0.9, 0.75]*diff(myXlims);
myYlims = get(gca,'Ylim');
ytr = myYlims(1)*ones(1,3) + [0.6, 0.6, 0.6+0.17*sqrt(3)]*diff(myYlims);
ctr = [0;1;2];
patch(xtr,ytr,ctr,'FaceVertexCData',ColorMatrix); % Make color patch
dxtr = 0.07;
dytr = 0.05;
trtxt(1) = text(xtr(1)-dxtr,ytr(1)-dytr,'Sand'); % Add labels to vertices
trtxt(2) = text(xtr(2)-0.5*dxtr,ytr(2)-dytr,'Silt');
trtxt(3) = text(xtr(3)-dxtr,ytr(3)+dytr,'Clay');
set(trtxt,'FontSize', 16)

hold off

% Set figure size and position
set(gcf, 'Units', 'Centimeters', 'Position', [10,10,24,14]);