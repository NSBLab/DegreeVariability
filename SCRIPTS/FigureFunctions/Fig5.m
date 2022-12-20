%% Fig 5a: Residual locations as a function of parcellation and pipeline

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
t1 = tiledlayout(2, length(parc), "TileSpacing", 'none', 'Padding', 'none', 'TileIndexing', 'columnmajor');

for jj = 1:length(parc)

    rois = roisAll{jj};
    currentSA = mean(PARC_SA{parc(jj)}, 1);

    current = DATA{ORDERED_INDS{parc(jj)}(pipe(jj)), THR_LEVELS};
    w = sum(current);
    r = calcResiduals(currentSA, w);

    nexttile(); axis off;
    brainplot(verts, faces, rois, r(1:end/2), 0, plasma);
    colorbar off;
    title({parc_name2{parc(jj)}, pipeline_titles{pipe(jj)}});

    nexttile(); axis off;
    brainplot(verts, faces, rois, r(1:end/2), 0, plasma);
    view([90 0]); colorbar off;


end

% linkaxes(ax(:));
% same_caxis(ax);
% title(t, currentDATA + "/" + parc_name{parc} + '/' + num2str(THR_LEVELS));
% print("v4_3a_" + currentDATA, '-dpng');

%% 5b: Strength and residuals Relationship

currentDATA = 'ADJGroupDen';
THR_LEVELS = 5;
parc = [5 3];
pipe = [7 5];


% ----

ax = gobjects(length(parc));

DATA = eval(currentDATA);

f = figure('Position', [149 125 550*length(parc) 500], 'DefaultAxesFontSize', 20);
t1 = tiledlayout(1, length(parc), "TileSpacing", 'none', 'Padding', 'loose');

SCALING = 10000;

for jj = 1:length(parc)

    currentSA = mean(PARC_SA{parc(jj)}, 1);

    current = DATA{ORDERED_INDS{parc(jj)}(pipe(jj)), THR_LEVELS};
    w = sum(current)/SCALING;
    r = calcResiduals(currentSA, w*SCALING);

    ax(jj) = nexttile();
    scatter(w, r, 100, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', .5); box on; axis square;

    bestfit();

    cp = corr(w', r');
    cs = corr(w', r','Type','Spearman');
    title("r = " + num2str(round(cp,4)), 'FontSize', 16, 'Units', 'normalized', 'Position', [0.50 0.9 0]);

    xtickformat('%.0f');
    ytickformat('%.1g');

   
%         ylabel(parc_name{parc(jj)});
%         ax(jj).YLabel.Visible = 'on';
    


end

xlabel(t1, 'Node Strength (\times10^4)', 'FontSize', 30); 
ylabel(ax(1), 'Residual Strength', 'FontSize', 30);
% linkaxes(ax(:));
% same_caxis(ax);
% title(t, currentDATA + "/" + parc_name{parc} + '/' + num2str(THR_LEVELS));
% print("v3_3b_" + currentDATA + "_" + num2str(THR_LEVELS), '-dpng');

%% Figure 5c: histograms of degree distribution for original and residuals

currentDATA = 'ADJGroupDen';
THR_LEVELS = 5;
parc = [5 3];
pipe = [7 5];

% -----

ax = gobjects(length(parc));

DATA = eval(currentDATA);

f = figure('Position', [149 125 550*length(parc) 500], 'DefaultAxesFontSize', 20);
t1 = tiledlayout(1, length(parc), "TileSpacing", 'compact', 'Padding', 'loose');

SCALING = 10000;

for jj = 1:length(parc)

    currentSA = mean(PARC_SA{parc(jj)}, 1);

    current = DATA{ORDERED_INDS{parc(jj)}(pipe(jj)), THR_LEVELS};
    w = sum(current)/SCALING;
    r = calcResiduals(currentSA, w*SCALING)/SCALING;

    ax(jj) = nexttile();
    
    hold on;
    

    [epdf, xi] = ksdensity(w);
    p1 = plot(xi, epdf, 'LineWidth', 2);
    
    [epdf, xi] = ksdensity(r);
    p2 = plot(xi, epdf, 'LineWidth', 2);

    h2 = histogram(r, 'Normalization', 'pdf', 'FaceColor', p2.Color, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);
    h1 = histogram(w, 'Normalization', 'pdf', 'FaceColor', p1.Color, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);

    yticks([]);
    box on;

%     xtickformat('%.0f');
%     ytickformat('%.1g');

   
%         ylabel(parc_name{parc(jj)});
%         ax(jj).YLabel.Visible = 'on';
    


end

legend('Original', 'Residuals');

xlabel(t1, 'Strength (\times10^4)', 'FontSize', 30); 
ylabel(ax(1), 'Relative Frequency', 'FontSize', 30);
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
    figure('Position', [ 100*ii, 80,1174, 729], 'DefaultAxesFontSize', 30)
    PlotMatrixWithLegend(corrsw(:, (10*ii-9):(10*ii)), ORDERED_MATRIX{parc(ii)}, LABELS{parc(ii)}, 'Correlation', threshs);
end


%%
% 
% figure('Position', [-1.3206e+03 221.8000 1.1808e+03 544], 'DefaultAxesFontSize', 18);
% tiledlayout(1, 2, 'TileSpacing','compact', 'Padding', 'compact')
% 
% nexttile();
% imagesc(corrsw);
% title({'Pearson Correlation', 'Strength vs Surface Area'});
% yticks([6, 17, 28, 39, 50]);
% yticklabels(currentDATALabels);
% 
% nexttile();
% imagesc(corrsd);
% title({'Pearson Correlation', 'Degree vs Surface Area'});
% yticks([6, 17, 28, 39, 50]); yticklabels([]); colorbar;
% 
% for ii = 1:2
%     nexttile(ii); hold on;
%     plot([10.5, 10.5], [0.5, 55.5], 'k');
%     plot([20.5, 20.5], [0.5, 55.5], 'k');
%     plot([30.5, 30.5], [0.5, 55.5], 'k');
%     plot([0.5, 40.5], [11.5, 11.5], 'k');
%     plot([0.5, 40.5], [22.5, 22.5], 'k');
%     plot([0.5, 40.5], [33.5, 33.5], 'k');
%     plot([0.5, 40.5], [44.5, 44.5], 'k');
%     set(gca,'YDir','normal');
% 
%     axis square;  caxis([0 1]);
%     xticks([5.5, 15.5, 25.5, 35.5]);
%     xticklabels(["82", "S220", "HCP", "S520"]);
% 
% end
% 
% nnz(corrsw(:)>0.8)/numel(corrsw)
% % figure; histogram(corrsw(:), 0:0.05:1);
% 
