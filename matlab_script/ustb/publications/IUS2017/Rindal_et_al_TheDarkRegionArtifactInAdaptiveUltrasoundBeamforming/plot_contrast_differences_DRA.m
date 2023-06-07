function [handle] = plot_contrast_differences_DRA(c,c_alt,c_alt_2,image,legend_1,legend_2,legend_3)
%PLOTLATERALLINE Plot lateral line from all images saved in image struct
%
C_all = [];

if nargin < 5
    legend_1 = 'Uncalibrated';
    legend_2 = 'Calibrated';
end

for i = 1:length(c) 
    C_all = [C_all; c(i) c_alt(i) c_alt_2(i)];
end
%%
    h101 = figure;
    subplot(211)
    bar(C_all)
    x_pos = [0.7 1.7 2.7 3.7 4.7 5.7 6.7 1 2 3 4 5 6 7 1.3 2.3 3.3 4.3 5.3 6.3 7.3];
    text(x_pos,double(C_all(:)),num2str(round(C_all(:)),'%d'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
    legend(legend_1,legend_2,legend_3);
       set(gca,'XTickLabel',image.tags)
       handle = gca;
end

