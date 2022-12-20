%% Fig 4a: Hub locations as a function of parcellation and pipeline

verts = lh_inflated_verts;
faces = lh_faces;
currentDATA = 'ADJGroupDen';
THR_LEVELS = 5;
parc = [5 3];
roisAll = { [lh_aparc], ...
    [Scha_parcs.lh_scha200], ...
    [lh_HCPMMP1], [Scha_parcs.lh_scha500]};

% ----

clear ax;
ax = gobjects(10, length(parc));

DATA = eval(currentDATA);

f = figure('Position', [500, 450, 1100, 160*length(parc) ], 'DefaultAxesFontSize', 30);
t1 = tiledlayout(length(parc), 10, "TileSpacing", 'none', 'Padding', 'none');

for jj = 1:length(parc)



    % re-index rois in right hemisphere
    rois = roisAll{jj};

    nodeLocations = zeros(max(rois), 3);
    for ii = 1:max(rois)
        nodeLocations(ii,:) = mean(verts(rois == ii, :), 1);
    end

%     currentSA = mean(PARC_SA{parc(jj)}, 1);

    for ii = 1:10


        current = DATA{ORDERED_INDS{parc(jj)}(ii), THR_LEVELS};
        w = sum(current);
        [~, worder] = sort(w, 'ascend');

        ax(ii,jj) = nexttile(t1); axis(ax(ii,jj), 'off');

        t3 = tiledlayout(t1, 2, 1, 'TileSpacing', 'none', 'Padding', 'none');
        t3.Layout.Tile = 10*(jj-1)+ii;
      

        nexttile(t3);
        brainplot(verts, faces, rois, w(1:end/2), 0, plasma);
        colorbar off; 

        nexttile(t3);
        brainplot(verts, faces, rois, w(1:end/2), 0, plasma);
        view([90 0]); colorbar off;

        

%         scatter(nodeLocations(worder,1), nodeLocations(worder,2), atan(currentSA(worder)/800)*50, ...
%             w(worder), 'filled', 'MarkerEdgeColor', 'black');
% 
%         axis equal; axis off;

        if ii == 1
            ylabel(ax(ii,jj), parc_name{parc(jj)});
            ax(ii,jj).YLabel.Visible = 'on';
        end

    end

end

% linkaxes(ax(:));
% same_caxis(ax);
% title(t, currentDATA + "/" + parc_name{parc} + '/' + num2str(THR_LEVELS));
scfw(1300);
% print("v4_3a_" + currentDATA, '-dpng');

%% 4b: Strength and SA Relationship

verts = [lh_verts; rh_verts];
currentDATA = 'ADJGroupDen';
THR_LEVELS = 5;
parc = [5 3];


% ----

clear ax;



ax = gobjects(10, 3);

DATA = eval(currentDATA);

f = figure('Position', [149 425 1.7376e+03 188*length(parc)], 'DefaultAxesFontSize', 30);
t1 = tiledlayout(length(parc), 10, "TileSpacing", 'none', 'Padding', 'loose');

for jj = 1:length(parc)

    currentSA = mean(PARC_SA{parc(jj)}, 1);

    for ii = 1:10


        current = DATA{ORDERED_INDS{parc(jj)}(ii), THR_LEVELS};
        w = sum(current);

        ax(ii,jj) = nexttile();
        scatter(currentSA, w, 10, 'filled'); box on; axis square;

        hold on;
        p = polyfit(currentSA, w,1);
        f = polyval(p,currentSA);
        plot(currentSA,f,'g', 'LineWidth', 2)

%         bestfit;
        bestfit('zeroIntercept', true);

        cp = corr(currentSA', w');
        cs = corr(currentSA', w','Type','Spearman');
        %         xlabel("Rho = " + num2str(round(cs,4)), 'FontSize', 12);

        xticklabels([]); xticks([]);
        yticklabels([]); yticks([]);

        if ii == 1
            ylabel(parc_name{parc(jj)});
            ax(ii,jj).YLabel.Visible = 'on';
        end

        if ii == 10
%             linkaxes(ax(1:10,jj));
        end

    end

end

xlabel(t1, 'Node Surface Area', 'FontSize', 30); ylabel(t1, 'Node Strength', 'FontSize', 30);
% linkaxes(ax(:));
% same_caxis(ax);
% title(t, currentDATA + "/" + parc_name{parc} + '/' + num2str(THR_LEVELS));
% print("v3_3b_" + currentDATA + "_" + num2str(THR_LEVELS), '-dpng');

%% Figure 4c: correlation coefficients for all parcellations, pipelines, metrics, thresholds

currentDATA = ["ADJGroupWei", "ADJGroupDen", "ADJGroupConDenAltered", "ADJGroupDst"];
currentDATALabels = ["Weight", "CV", "Consistency", "Distance Bins"];
THR_LEVELS = 1:11;
parc = [1 5 3 6];



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

figure('Position', [-1.3206e+03 221.8000 1.1808e+03 544], 'DefaultAxesFontSize', 18);
tiledlayout(1, 2, 'TileSpacing','compact', 'Padding', 'compact')

nexttile();
imagesc(corrsw);
title({'Pearson Correlation', 'Strength vs Surface Area'});
yticks([6, 17, 28, 39, 50]);
yticklabels(currentDATALabels);

nexttile();
imagesc(corrsd);
title({'Pearson Correlation', 'Degree vs Surface Area'});
yticks([6, 17, 28, 39, 50]); yticklabels([]); colorbar;

for ii = 1:2
    nexttile(ii); hold on;
    plot([10.5, 10.5], [0.5, 55.5], 'k');
    plot([20.5, 20.5], [0.5, 55.5], 'k');
    plot([30.5, 30.5], [0.5, 55.5], 'k');
    plot([0.5, 40.5], [11.5, 11.5], 'k');
    plot([0.5, 40.5], [22.5, 22.5], 'k');
    plot([0.5, 40.5], [33.5, 33.5], 'k');
    plot([0.5, 40.5], [44.5, 44.5], 'k');
    set(gca,'YDir','normal');

    axis square;  caxis([0 1]);
    xticks([5.5, 15.5, 25.5, 35.5]);
    xticklabels(["82", "S220", "HCP", "S520"]);

end

nnz(corrsw(:)>0.8)/numel(corrsw)
% figure; histogram(corrsw(:), 0:0.05:1);

