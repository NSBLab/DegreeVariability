%% Fig 4a: Hub locations as a function of parcellation and pipeline

verts = lh_inflated_verts;
faces = lh_faces;
currentDATA = 'ADJGroupDen';
THR_LEVELS = 5;
parc = [5 3];
pipe = [7 5];
roisAll = {  [Scha_parcs.lh_scha200], [lh_HCPMMP1] };

% ----

clear ax;
ax = gobjects(length(parc));

DATA = eval(currentDATA);

f = figure('Position', [500, 450, 480*length(parc), 700 ], 'DefaultAxesFontSize', 30);
% t1 = tiledlayout(2, length(parc), "TileSpacing", 'none', 'Padding', 'none', 'TileIndexing', 'columnmajor');
t1 = tiledlayout(2, length(parc), "TileSpacing", 'none', 'Padding', 'none');

for jj = 1:length(parc)

    rois = roisAll{jj};


    current = DATA{ORDERED_INDS{parc(jj)}(pipe(jj)), THR_LEVELS};
    w = sum(current);


    nexttile(); axis off;
    brainplot(verts, faces, rois, w(1:end/2), 0, plasma);
    colorbar off;
%     title({parc_name2{parc(jj)}, pipeline_titles{pipe(jj)}});

    nexttile(); axis off;
    brainplot(verts, faces, rois, w(1:end/2), 0, plasma);
    view([90 0]); colorbar off;


end

% linkaxes(ax(:));
% same_caxis(ax);
% title(t, currentDATA + "/" + parc_name{parc} + '/' + num2str(THR_LEVELS));
% print("v4_3a_" + currentDATA, '-dpng');

%% 4b: Strength and SA Relationship

% currentDATA = 'ADJGroupDen';
% THR_LEVELS = 5;
% parc = [5 3];
% pipe = [7 5];


% ----

ax = gobjects(length(parc));

DATA = eval(currentDATA);

% f = figure('Position', [149 125 550*length(parc) 500], 'DefaultAxesFontSize', 20);
% t1 = tiledlayout(1, length(parc), "TileSpacing", 'none', 'Padding', 'loose');

f = figure('Position', [149 125 480 700], 'DefaultAxesFontSize', 24);
t1 = tiledlayout( length(parc), 1, "TileSpacing", 'tight', 'Padding', 'tight');

SCALING = 10000;

for jj = 1:length(parc)

    currentSA = mean(PARC_SA{parc(jj)}, 1);

    current = DATA{ORDERED_INDS{parc(jj)}(pipe(jj)), THR_LEVELS};
    w = sum(current)/SCALING;

    ax(jj) = nexttile();
    scatter(currentSA, w, 100, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5); box on; axis square;

    bestfit('plotOptions', {'r', 'LineWidth', 4});
%     hold on;
%     p = polyfit(currentSA, w,1);
%     f = polyval(p,currentSA);
%     plot(currentSA,f,'r', 'LineWidth', 2)

    cp = corr(currentSA', w');
    cs = corr(currentSA', w', 'Type','Spearman');
    title("r = " + num2str(round(cp,4)), 'FontSize', 24, 'Units', 'normalized', 'Position', [0.50 0.85 0]);

    temp = xticks; xticks([temp(1), temp(end)]);
    temp = yticks; yticks([temp(1), temp(end)]);

   
%         ylabel(parc_name{parc(jj)});
%         ax(jj).YLabel.Visible = 'on';
    


end

xlabel(t1, 'Node Surface Area', 'FontSize', 30); 
ylabel(t1, 'Node Strength (\times10^4)', 'FontSize', 30);
% linkaxes(ax(:));
% same_caxis(ax);
% title(t, currentDATA + "/" + parc_name{parc} + '/' + num2str(THR_LEVELS));
% print("v3_3b_" + currentDATA + "_" + num2str(THR_LEVELS), '-dpng');

%% Figure 4c: correlation coefficients for example parcellations, pipelines,  thresholds

currentDATA = ["ADJGroupDen"];
THR_LEVELS = 1:11;
parc = [5 3];

% ----

corrsw = zeros(length(currentDATA)*length(THR_LEVELS), length(parc)*10);
corrsd = zeros(size(corrsw));

for parci = 1:length(parc)

    currentSA = mean(PARC_SA{parc(parci)}, 1);

    for metri = 1:length(currentDATA)

        DATA = eval(currentDATA(metri));

        for pipei = 1:10

            for thri = 1:length(THR_LEVELS)

                current = DATA{ORDERED_INDS{parc(parci)}(pipei), thri};
                w = sum(current);
                d = sum(~~current);

                cpw = corr(currentSA', w');
                cpd = corr(currentSA', d');


                corrsw(...
                    (metri-1)*length(THR_LEVELS) + thri, ...
                    (parci-1)*10 + pipei                    ) = cpw;

                corrsd(...
                    (metri-1)*length(THR_LEVELS) + thri, ...
                    (parci-1)*10 + pipei                    ) = cpd;

            end
        end
    end
end


%%

% figure('Position', [-1.3206e+03 221.8000 1.1808e+03 544], 'DefaultAxesFontSize', 18);
% tiledlayout(1, length(parc), 'TileSpacing','compact', 'Padding', 'compact')

for ii = 1:length(parc)
    figure('Position', [ 100*ii, 80, 700, 420], 'DefaultAxesFontSize', 24)
    PlotMatrixWithLegend(corrsw(:, (10*ii-9):(10*ii)), ORDERED_MATRIX{parc(ii)}, LABELS{parc(ii)}, 'Correlation', threshs, 'caxis', [0 1], 'cbarOptions', {'XTick', 0:0.25:1}, 'cbarPosition', [.15 .35 .05 .6]);
end

