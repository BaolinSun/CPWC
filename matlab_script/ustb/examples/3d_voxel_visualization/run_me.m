% example of 3d US data visualization uing vol3d

clear all
close all
clc

load('3d_data.mat')

dyn  = 45;
gain = -100;

use_painters = false; % set to true for best visualization effect, but really slow           

% === scan conversion  ===
az = linspace(bmode.geo.axes(1), bmode.geo.axes(4), bmode.geo.gridSize(2));
el = linspace(bmode.geo.axes(2), bmode.geo.axes(5), bmode.geo.gridSize(3));
r  = linspace(bmode.geo.axes(3), bmode.geo.axes(6), bmode.geo.gridSize(1));

zq = r;
xq = linspace(r(end)*sin(az(1)), r(end)*sin(az(end)), 100);
yq = linspace(r(end)*sin(el(1)), r(end)*sin(el(end)), 100);

[R, AZ, EL]    = ndgrid(r, az, el);
[Xq, Yq, Zq]   = ndgrid(xq, yq, zq);
[AZq, ELq, Rq] = cart2sph(Zq, Xq, Yq);


data = bmode.data(:,:,:, 1);
data = interpn(R, AZ, EL, data, Rq, AZq, ELq, 'linear');

% === compression ===
data = 20*log10(data);

% === the intensity is converted into transparency values instead of grascale values ===

ColorInd = reshape(ind2rgb(reshape(cast(256*(data-gain)/dyn, 'uint8'), [numel(data),1]), gray(256)), [size(data), 3]);

AlphaInd = ones(size(data));
AlphaInd(data<gain)   = 0;
AlphaInd(isnan(data)) = 0;


h = vol3d('CData', ColorInd, 'Alpha', AlphaInd, ...
    'XData', [min(xq), max(xq)], 'YData', [min(yq), max(yq)], 'ZData', [min(zq), max(zq)]);
grid on
set(gcf, 'Color', 'k')
set(gca, 'Color', 'none')
set(gca, 'XColor', 'w')
set(gca, 'YColor', 'w')
set(gca, 'ZColor', 'w')
set(gca, 'box', 'on')
set(gca, 'GridLineStyle', ':')
set(gca, 'GridColor', 'w')
set(gca, 'GridAlpha', 0.75)
set(gca, 'Xdir', 'Reverse')
set(gca, 'Zdir', 'Reverse')
axis equal tight

if use_painters
    set(gcf, 'renderer', 'painters')
else
    set(gcf, 'renderer', 'openGL')
end

