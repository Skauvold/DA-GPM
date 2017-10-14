% A script for creating Fig. 1

clear; close all; clc

% Load data
load('./data/ANSCaseEnKF_output.mat');

%%

% Grid size
nx = 110;
ny = 87;

% Time steps
steps = 50;
tlen = steps+1;

% Organize state variables
final_analysis = Xa;

ZS_final = final_analysis(2*tlen+1:end,:);
ZS_array = reshape(mean(ZS_final,2),nx,ny,[]);

Z_final = ZS_array(:,:,1:steps+1);
S1_final = ZS_array(:,:,steps+2:2*steps+1);
S2_final = ZS_array(:,:,2*steps+2:3*steps+1);
S3_final = ZS_array(:,:,3*steps+2:4*steps+1);

zOut = Z_final;
s1Out = S1_final;
s2Out = S1_final;
s3Out = S1_final;

%%

% Well location on 2D grid
x_Tunalik1 = 69-3;
y_Tunalik1 = 67-2;

length_x = 0.25;
length_y = 0.25;

f_Tunalik1 = 0.55;

% Number of interpolation points
nq = 50;

% Set up interpolation axis
slicex_start = x_Tunalik1/nx - f_Tunalik1*length_x;
slicex_end = x_Tunalik1/nx + (1-f_Tunalik1)*length_x;
slicey_start = y_Tunalik1/ny - f_Tunalik1*length_y;
slicey_end = y_Tunalik1/ny + (1-f_Tunalik1)*length_y;

xq = linspace(slicex_start,slicex_end,nq);
yq = linspace(slicey_start,slicey_end,nq);
[Xv,Yv] = meshgrid(linspace(0,1,nx),linspace(0,1,ny));

zq = [];
Tp = [0.5,0.5,0,0;    0,0,0.6,0.4;    0,0,0.2,0.8];

% Interpolate state variables onto slice
for i = 1:51
    zq(i,:) = interp2(Xv,Yv,zOut(:,:,i)',xq,yq);
    if i <= 50
        s1qi = interp2(Xv,Yv,s1Out(:,:,i)',xq,yq);
        s2qi = interp2(Xv,Yv,s2Out(:,:,i)',xq,yq);
        s3qi = interp2(Xv,Yv,s3Out(:,:,i)',xq,yq);
        pqi = eta2p([s1qi;s2qi;s3qi]);
        pqi = [pqi;ones(1,size(pqi,2))-sum(pqi,1)];
        pqi = Tp*pqi;
        pqi = pqi./repmat(sum(pqi,1),3,1);
        p1q(i,:) = pqi(1,:);
        p2q(i,:) = pqi(2,:);
        p3q(i,:) = pqi(3,:);
    end
end

% Rescale along-slice axis
Lx = 15*nx; % Domain size [km]
Ly = 15*ny;
slice_length = sqrt(Lx*Lx*(slicex_end-slicex_start)^2 + ...
    Ly*Ly*(slicey_end-slicey_start)^2);
slice_scale = linspace(0,slice_length,nq);

% Set colours
Color1 = [255, 0, 0]./255; % red
Color2 = [255, 246, 0]./255; % yellow
Color3 = [60, 206, 30]./255; % green
ColorMat = [Color1; Color2; Color3];

% Make figure
figure()
hold on
for layer = 1:50
    % Get bottom and top curves
    zbot = zq(layer,:);
    ztop = zq(layer+1,:);
    for iCell = 1:(nq-1)
        % Compute color
        plc = [p1q(layer,iCell);p2q(layer,iCell);p3q(layer,iCell)];
        colorlc = plc'*ColorMat;
        % Prepare cell edges
        zlc_bot = zbot(iCell:(iCell+1));
        zlc_top = ztop(iCell:(iCell+1));
        zlc = [zlc_bot,fliplr(zlc_top)];
        xlc_bot = slice_scale(iCell:(iCell+1));
        xlc = [xlc_bot,fliplr(xlc_bot)];
        % Plot patch of cell
        patch(xlc,zlc,colorlc, 'EdgeColor', 'None');
    end
end
plot(slice_scale,zq(1,:),'k-');
plot(slice_scale,zq(end,:),'k-');

xlim(slice_scale([1,end]))
ylim([-3200, 2000])

% Add colour legend to top right part of figure
myXlims = get(gca,'Xlim');
xtr = myXlims(1)*ones(1,3) + [0.6, 0.9, 0.75]*diff(myXlims);
myYlims = get(gca,'Ylim');
ytr = myYlims(1)*ones(1,3) + [0.6, 0.6, 0.6+0.17*sqrt(3)]*diff(myYlims);
ctr = [0;1;2];
dytr = 250;
dxtr = 20;
patch(xtr,ytr,ctr,'FaceVertexCData',ColorMat);
trtxt(1) = text(xtr(1)-dxtr,ytr(1)-dytr,'Sand');
trtxt(2) = text(xtr(2)-dxtr,ytr(2)-dytr,'Silt');
trtxt(3) = text(xtr(3)-dxtr,ytr(3)+dytr,'Clay');
set(trtxt,'FontSize',16)
set(gca,'FontSize',14)
xlabel('Distance along section [km]')
ylabel('Elevation [m]')

% Add direction labels to top edge of figure
text(myXlims(1), myYlims(2)+0.05*diff(myYlims),'SW', 'FontSize', 16)
text(myXlims(2)-0.05*diff(myXlims),myYlims(2)+0.05*diff(myYlims),'NE', 'FOntSize', 16)

box on
plot(f_Tunalik1*slice_scale(end)*ones(1,2), myYlims, 'k:')

% Add "mini-map" to bottom left part of figure
minimaphandle = axes('Position',[.17 .17 .22 .22]);
hold on
contour(Xv,Yv,zOut(:,:,51)', [-3000:500:3000], 'k-')
plot(xq, yq, 'k-', 'LineWidth',1.6)
plot(x_Tunalik1/nx, y_Tunalik1/ny, 'ko', 'LineWidth', 1.6)
annotation('arrow',(0.175+0.13*0.22)*[1,1],[0.17+0.4*0.22, 0.17+0.9*0.22])
text(0.1,0.6,'N', 'FontSize', 16)


hold off
view([0,90])
set(minimaphandle, 'XTick',[], 'YTick', [])

hold off
box on