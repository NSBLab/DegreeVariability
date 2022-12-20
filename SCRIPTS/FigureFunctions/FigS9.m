%% Heatmaps
% for a given parc, over all pipelines




groupName = {'ADJGroupWei', 'ADJGroupDen', 'ADJGroupConDenAltered'};
parc = 5;
pipeline = 1:10;
thrLevels = [1 3 5 7 9 11];


% f = figure;%("Position",[1.5378e+03 940.2 2.0464e+03 600]);
%t = tiledlayout(2, 5, 'TileSpacing', 'Tight');
ax = [];
heatmapData.group = [];
out = [];


figure('DefaultAxesFontSize', 5);
tl = tiledlayout(5, 2, 'TileIndexing', 'columnmajor', 'TileSpacing', 'tight');



for jj = 1:length(pipeline)

    temp = [1, length(thrLevels), length(groupName)];
    [a, b, c] = ind2sub(temp, 1:prod(temp));

    for ii = 1:prod(temp)

        currPipelineNo = ORDERED_INDS{parc}(pipeline(jj));

        heatmapData.group = eval(groupName{c(ii)});

        current = heatmapData.group{currPipelineNo, thrLevels(b(ii))};

        out(ii, :) = sum(current);

    end

    ax=nexttile();

    CSTR = partialcorr(out', mean(PARC_SA{parc})', 'Type','Spearman');
    imagesc(CSTR); axis square

    hold on;
    xline(6.5); xline(12.5); yline(6.5); yline(12.5);

    % ticks v1
%     xticks([2,5,8,11,14,17]); xticklabels(thr_strings_density(repmat(1:10:11, 1, 3))+"%");
%     xtickangle(45); xlabel('Density');
%     yticks([3.5, 9.5, 15.5]); yticklabels({'Weight', 'CV', 'Consistency'});

    % ticks v2
    yticks([2,3.5,5,8,9.5,11,14,15.5,17]); yticklabels({'5%', 'Weight          ', '30%', '5%', 'CV          ', '30%', '5%', 'Consistency          ', '30%'});
    set(gca, 'ticklength', [0 0]); xticks([]);

    % ticks v3
%     xticks([2,5,8,11,14,17]); xticklabels(thr_strings_density(repmat(1:10:11, 1, 3)));
%     xlabel(tl, 'Density (%)');
%     yticks([3.5, 9.5, 15.5]); yticklabels({'Weight', 'CV', 'Consistency'});
%     set(gca, 'ticklength', [0 0]);

    title(pipeline_titles(jj));
    %     set(gca, 'FontSize', 5);

    caxis([0 1])
    colormap(flipud(make_cmap('orangered',250,30,0)));

    min(CSTR(:))

end
scfw(600); scfh(850);


%%

matrixLabels = [repmat(ORDERED_MATRIX{parc}, 1, length(groupName)*length(thrLevels)); ...
    reshape(repmat(0:(length(groupName)-1),length(pipeline)*length(thrLevels),1),1, []) ];

labelInfo = [LABELS{parc};
    {'ThrMetric:{\color[rgb]{0.729412,0.729412,0.729412}Weight}/{\color[rgb]{0.250980,0.250980,0.250980}CV}/{\color[rgb]{0.956863,0.647059,0.509804}Con}'}];

% CSTR = corr(out', 'Type','Spearman');
CSTR = partialcorr(out', mean(PARC_SA{parc})', 'Type','Spearman');
RunClusterPipelineProp(CSTR,2,1,matrixLabels,labelInfo, ...
    [parc_name{parc}, ' rank correlation']);
set(gcf, 'Position', [-1535 107.4000 1536 740.8000]);