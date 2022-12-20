%% Heatmaps
% for a given parc, over all pipelines




groupName = {'ADJGroupWei', 'ADJGroupDen', 'ADJGroupConDenAltered'};
parc = 5;
pipeline = 1:10;
thrLevels = 7;





% f = figure;%("Position",[1.5378e+03 940.2 2.0464e+03 600]);
%t = tiledlayout(2, 5, 'TileSpacing', 'Tight');
ax = [];
heatmapData.group = [];
out = [];

temp = [length(pipeline), length(thrLevels), length(groupName)];
[a, b, c] = ind2sub(temp, 1:prod(temp));

for ii = 1:prod(temp)

    currPipelineNo = ORDERED_INDS{parc}(pipeline(a(ii)));

    heatmapData.group = eval(groupName{c(ii)});

    current = heatmapData.group{currPipelineNo, thrLevels(b(ii))};

    out(ii, :) = sum(current);
%     out(ii, :) = sum(logical(current));
end

%%

matrixLabels = [repmat(ORDERED_MATRIX{parc}, 1, length(groupName)*length(thrLevels)); ...
    reshape(repmat(0:(length(groupName)-1),length(pipeline)*length(thrLevels),1),1, []) ];

labelInfo = [LABELS{parc};
    {'ThrMetric:{\color[rgb]{0.729412,0.729412,0.729412}Weight}/{\color[rgb]{0.250980,0.250980,0.250980}CV}/{\color[rgb]{0.956863,0.647059,0.509804}Con}'}];

% CSTR = corr(out', 'Type','Spearman');
CSTR = partialcorr(out', mean(PARC_SA{parc})', 'Type','Spearman');
[~,clusterAllocations] = RunClusterPipelineProp(CSTR,2,1,matrixLabels,labelInfo, ...
    [parc_name2{parc}, ' rank correlation']);
set(gcf, 'Position', [1000 0100 1000 0900]);

%%

figure('DefaultAxesFontSize', 18); 
n = max(clusterAllocations);
tl = tiledlayout(n+1, 1, 'TileSpacing', 'tight', 'Padding', 'none');
nexttile();
binEdges = linspace(-0.5, 1, 16);
histogram(CSTR(~eye(size(CSTR))), binEdges);
box on;
xlim([-0.55 1.05]);
xticklabels([]);
title('All Correlations')

for ii = 1:n
    nexttile();
    temp = CSTR(clusterAllocations==ii, clusterAllocations==ii);
    histogram(temp(~eye(size(temp))), binEdges);
    box on;
    xlim([-0.55 1.05]);
    xticklabels([]);
    title(sprintf('Cluster %i', ii));
end

xticklabels([-.5, 0, .5, 1]); xtickangle(45)
xlabel(tl, '\rho', 'FontSize', 18);
ylabel(tl, 'Count', 'FontSize', 18);
scfw(250); scfh(500);

%%

% print("v3_5_" + parc_name{parc} + "_" + num2str(thrLevels), '-dpng', '-r600');

% get location of dissimilar pipelines
[minrow, mincol] = find(CSTR == min(CSTR(:)), 1);

% get pipeline details
f = @(x) ([a(x), b(x), c(x)]);
locs(1,:) = f(minrow);
locs(2,:) = f(mincol);

figure;
layout1 = tiledlayout(1, 2, 'TileSpacing', 'tight', 'Padding', 'none');
for ii = 1:size(locs, 1)

    layout2 = tiledlayout(layout1, 2, 2, 'TileSpacing', 'none', 'Padding', 'none', 'TileIndexing', 'columnmajor');
    layout2.Layout.Tile = ii;

    temp = eval(groupName{c(locs(ii,3))});
    currentData = temp{ORDERED_INDS{parc}(locs(ii,1)), locs(ii,2)};
    wrank = tiedrank(sum(currentData));

    for jj = 1:4
        if jj < 2.5
            verts = lh_verts; faces = lh_faces;
            rois = Scha_parcs.lh_scha200;
            w = wrank(1:100);
        else
            verts = rh_verts; faces = rh_faces;
            rois = Scha_parcs.rh_scha200;
            w = wrank(101:200);
        end

        nexttile(layout2);
        
        brainplot(verts, faces, rois, w, 0, plasma); colorbar off;

        if jj > 1.5 && jj < 3.5
            set(gca, 'View', [90 0]);
        end

    end

end

scfw(1500); scfh(500);