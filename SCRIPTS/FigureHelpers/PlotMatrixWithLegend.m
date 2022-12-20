function PlotMatrixWithLegend(DATA,PROCESSING_MATRIX,PROCESSING_MATRIX_LABELS,colorbar_label,yticklabel_data,varargin)

% Plots a matrix of data for all different processing pipelines

% ---

ip = inputParser;
addRequired(ip, 'DATA');
addRequired(ip, 'PROCESSING_MATRIX');
addRequired(ip, 'PROCESSING_MATRIX_LABELS');
addRequired(ip, 'colorbar_label');
addRequired(ip, 'yticklabel_data');

addOptional(ip, 'ylabel_name', 'Threshold');

addParameter(ip, 'cbarStatus', 'on');
addParameter(ip, 'cbarOptions', {});
addParameter(ip, 'cmap', viridis);
addParameter(ip, 'caxis', "auto");
addParameter(ip, 'cbarPosition', [.125 .35 .05 .6]);

addParameter(ip, 'plotParent', gcf);
addParameter(ip, 'ylabelStatus', 'on');

parse( ip , DATA,PROCESSING_MATRIX,PROCESSING_MATRIX_LABELS,colorbar_label,yticklabel_data, varargin{:} );

% ---

PanelAPos = [.5 .35 .5 .6];
PanelBPos = [.5 .10 .5 .25];

n = size(PROCESSING_MATRIX,2)+1;

ax1 = axes(ip.Results.plotParent, 'Position',PanelAPos);

imagesc(DATA)

yticks(1:2:11); 
xticks(1:11); xticklabels([]);

if (strcmp(ip.Results.ylabelStatus, 'on'))
    yticklabels(yticklabel_data(1:2:end)*100 + "%")
    ylabel(ip.Results.ylabel_name)
else
    yticklabels([]);
end

colormap(ax1, ip.Results.cmap); 
caxis(ip.Results.caxis);

if (strcmp(ip.Results.cbarStatus, 'on'))
    c = colorbar(ax1,'Position',ip.Results.cbarPosition, ip.Results.cbarOptions{:});
    c.Label.String = colorbar_label;
end



ax2 = axes(ip.Results.plotParent, 'Position',PanelBPos);

Color1 = [186,186,186]./255;
Color2 = [64,64,64]./255;
Color3 = [244,165,130]./255;
Color4 = [171,217,233]./255;
Color5 = [69,117,180]./255;
Color6 = [49,54,149]./255;

imagesc(PROCESSING_MATRIX)
yticks(1:length(PROCESSING_MATRIX_LABELS));
ytickangle(0)
xticks(1:n-1)
xtickangle(90);
hold on

if (strcmp(ip.Results.ylabelStatus, 'on'))
    yticklabels(PROCESSING_MATRIX_LABELS);
    set(gca, 'YTickLabel', PROCESSING_MATRIX_LABELS);
else
    yticklabels([]);
end

for i = 0:1:n
    plot([0 n],[i-.5 i-.5],'k')
end

for i = 0:n
    plot([i-.5 i-.5],[0 n],'k')
end

xlabel('Pipeline')

%colormap(ax2,[Color1; Color2; Color3; Color4; Color5; Color6])
% colormap(ax2,[Color1; Color2; Color3; Color4; Color5])
colormap(ax2,[Color1; Color2; Color3])