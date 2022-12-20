%% Heatmaps
% for a given parc, over all pipelines




groupName = {'ADJGroupWei', 'ADJGroupDen', 'ADJGroupConDenAltered'};
parc = [1 5 3 6];
pipeline = 1:10;
thrLevels = 7;
%-----------------------


figure;
tl1 = tiledlayout(length(parc), 1);

for kk = 1:length(parc)
    % f = figure;%("Position",[1.5378e+03 940.2 2.0464e+03 600]);
    %t = tiledlayout(2, 5, 'TileSpacing', 'Tight');
    ax = [];
    heatmapData.group = [];
    out = [];

    temp = [length(pipeline), length(thrLevels), length(groupName)];
    [a, b, c] = ind2sub(temp, 1:prod(temp));

    for ii = 1:prod(temp)

        currPipelineNo = ORDERED_INDS{parc(kk)}(pipeline(a(ii)));

        heatmapData.group = eval(groupName{c(ii)});

        current = heatmapData.group{currPipelineNo, thrLevels(b(ii))};

        out(ii, :) = sum(current);
        %     out(ii, :) = sum(logical(current));
    end

    matrixLabels = [repmat(ORDERED_MATRIX{parc(kk)}, 1, length(groupName)*length(thrLevels)); ...
        reshape(repmat(0:(length(groupName)-1),length(pipeline)*length(thrLevels),1),1, []) ];

    % labelInfo = [LABELS{parc(kk)};
    %     {'ThrMetric:{\color[rgb]{0.729412,0.729412,0.729412}Weight}/{\color[rgb]{0.250980,0.250980,0.250980}CV}/{\color[rgb]{0.956863,0.647059,0.509804}Con}'}];

    CSTR = partialcorr(out', mean(PARC_SA{parc(kk)})', 'Type','Spearman');
    % [~,clusterAllocations] = RunClusterPipelineProp(CSTR,2,0,matrixLabels,labelInfo, ...
    %     [parc_name2{parc}, ' rank correlation']);

    tract = matrixLabels(3, :);
    CSTR(~~tril(ones(size(CSTR)))) = nan;


    % plotting
    tl2 = tiledlayout(tl1, 1, 4);
    tl2.Layout.Tile = kk;

    for jj = 1:4

        switch jj
            case 1; out = CSTR; t = 'All pipelines'; % Overall histogram
            case 2; out = (CSTR(tract==0, tract==0)); t = 'Deterministic only'; % within det
            case 3; out = (CSTR(tract==1, tract==1)); t = 'Probabilistic only';% within prob
            case 4; out = (CSTR(tract ~= tract')); % implicit expansion
                t = {'Deterministic vs', 'Probabilistic'}; % det vs prob
        end

        ax(end+1) = nexttile(tl2);

%         h = histogram(out(:), linspace(-1, 1, 101), 'Normalization', 'pdf', 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
        h = histogram(out(:), linspace(-1, 1, 101), 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
        xlim([-1 1]); yticks([]); xticks([-1:0.5:1]); xticklabels([]);

        [ks, xi] = ksdensity(out(:), h.BinEdges, 'Boundary', 'reflection', 'Function', 'pdf');
        % plot ksdensity so that the area under the histogram and the
        % density plot are the same
        hold on; plot(xi, ks*sum(h.Values)/sum(ks), 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);
        set(gca, 'TickLength', [.05 0]);

        switch kk
            case 1; title(t);
            case length(parc); xticklabels({-1, [], 0,[] , 1});
        end

    end

    ylabel(tl2, parc_name2(parc(kk)));

linkaxes(ax, 'y');

end
xlabel(tl1, '\rho');    
ylabel(tl1, 'Relative Frequency');