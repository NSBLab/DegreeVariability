%% Figure 3a: skewness and kurtosis for strength distributions for all parcellations, pipelines, metrics, thresholds

currentDATA = ["ADJGroupWei", "ADJGroupDen", "ADJGroupConDenAltered", "ADJGroupDst"];
currentDATALabels = ["Weight", "CV", "Consistency", "Distance Bins"];
THR_LEVELS = 1:11;
parc = [1 5 3 6];



% ----

skeww = zeros(length(currentDATA)*length(THR_LEVELS), length(parc)*10);
skewwr = zeros(size(skeww));
skewd = zeros(size(skeww));
kurtw = zeros(size(skeww));
kurtwr = zeros(size(skeww));
kurtd = zeros(size(skeww));

for parci = 1:length(parc)

    currentSA = mean(PARC_SA{parc(parci)}, 1);

    for metri = 1:length(currentDATA)

        DATA = eval(currentDATA(metri));

        for pipei = 1:10

            for thri = 1:length(THR_LEVELS)

                current = DATA{ORDERED_INDS{parc(parci)}(pipei), thri};
                w = sum(current);
                d = sum(logical(current));

                p = polyfit(currentSA, w,1);
                wHat = polyval(p,currentSA);
                wr = w-wHat; % weight_residuals

                currentX = (metri-1)*length(THR_LEVELS) + thri;
                currentY = (parci-1)*10 + pipei;

                skeww(currentX, currentY) = skewness(w);
                skewwr(currentX, currentY) = skewness(wr); 

                skewd(currentX, currentY) = skewness(d);

                kurtw(currentX, currentY) = kurtosis(w);
                kurtwr(currentX, currentY) = kurtosis(wr);

                kurtd(currentX, currentY) = kurtosis(d);

            end
        end
    end
end

%%

figure('Position', [-1.3206e+03 221.8000 1.1808e+03 544], 'DefaultAxesFontSize', 18);
tiledlayout(2, 2, 'TileSpacing','compact', 'Padding', 'compact')

nexttile();
imagesc(skeww);
title({'Skewness', 'Strength'});
yticks([6, 17, 28, 39, 50]); caxis([0 2]);
yticklabels(currentDATALabels);

nexttile();
imagesc(skewwr);
title({'Skewness', 'Strength Residuals'});
yticks([6, 17, 28, 39, 50]); caxis([0 2]);
yticklabels([]); colorbar;

nexttile();
imagesc(kurtw-3);
title({'Excess Kurtosis', 'Strength'});
yticks([6, 17, 28, 39, 50]); caxis([0 6]);
yticklabels(currentDATALabels);

nexttile();
imagesc(kurtwr-3);
title({'Excess Kurtosis', 'Strength Residuals'});
yticks([6, 17, 28, 39, 50]); caxis([0 6]);
yticklabels([]); colorbar;

for ii = 1:4
    nexttile(ii); hold on;
    plot([10.5, 10.5], [0.5, 55.5], 'k');
    plot([20.5, 20.5], [0.5, 55.5], 'k');
    plot([30.5, 30.5], [0.5, 55.5], 'k');
    plot([0.5, 40.5], [11.5, 11.5], 'k');
    plot([0.5, 40.5], [22.5, 22.5], 'k');
    plot([0.5, 40.5], [33.5, 33.5], 'k');
    plot([0.5, 40.5], [44.5, 44.5], 'k');
    set(gca,'YDir','normal');

%     axis square;  %caxis([0 1]);
    xticks([5.5, 15.5, 25.5, 35.5]);
    xticklabels(["82", "S220", "HCP", "S520"]);

end

scfh(1000); scfw(1000);
% print('v3_4a', '-dpng', '-r600');



%%
currentDATA = ["ADJGroupWei", "ADJGroupDen", "ADJGroupConDen", "ADJGroupCon", "ADJGroupDst"];
currentDATALabels = ["Weight", "CV", "Consistency+Density", "Consistency", "Distance Bins"];
THR_LEVELS = 1:11;
parc = [1 5 3 6];



for parci = 1:length(parc)

    currentSA = mean(PARC_SA{parc(parci)}, 1);

    for metri = 1:length(currentDATA)

        DATA = eval(currentDATA(metri));

        for pipei = 1:10

            for thri = 1:length(THR_LEVELS)


                current = DATA{ORDERED_INDS{parc(parci)}(pipei), thri};
                w = sum(current);


                [~,~,timr(... % tail index mild right tail
                    (metri-1)*length(THR_LEVELS) + thri, ...
                    (parci-1)*10 + pipei                    ),...
                    tier(... % tail index extreme right tail
                    (metri-1)*length(THR_LEVELS) + thri, ...
                    (parci-1)*10 + pipei                    )] = getFenceIndices(w);

            end
        end
    end
end


%%


figure('Position', [-1.3206e+03 221.8000 1.1808e+03 544], 'DefaultAxesFontSize', 18);
tiledlayout(1, 2, 'TileSpacing','compact', 'Padding', 'compact')

nexttile();
imagesc(timr);
title({'Intermediate Right Tail', 'Strength'});
yticks([6, 17, 28, 39, 50]);
yticklabels(currentDATALabels); colorbar; caxis([0 0.0339]);

nexttile();
imagesc(tier);
title({'Extreme Right Tail', 'Strength'});
yticks([6, 17, 28, 39, 50]); yticklabels([]);colorbar;

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

    axis square;
    xticks([5.5, 15.5, 25.5, 35.5]);
    xticklabels(["82", "S220", "HCP", "S520"]);

end


scfh(500); scfw(1000);
print('v3_4b', '-dpng', '-r600');

%% 3c example of kurtosis and skewness

parc = 3;
pipeline = 6;
currentDATA = "ADJGroupConDen";

figure;
tiledlayout(1, 11);
DATA = eval(currentDATA);

clear ax;
for ii = 1:11

    ax(ii) = nexttile();


    current = DATA{ORDERED_INDS{parc}(pipeline), THR_LEVELS(ii)};
    [epdf, xi] = ksdensity(sum(current~=0), 'Function', 'pdf');
    plot(xi, epdf, 'LineWidth', 2);


end

linkaxes(ax);