parci = 5;
pipei = 7;

DATA = ADJGroupConDenAltered;
thri = 7;

surfaceAreas = mean(PARC_SA{parci}, 1); whos surfaceAreas
edgeWeights = DATA{ORDERED_INDS{parci}(pipei), thri}; whos edgeWeights

%%
figure('Position', [400 050 500 1000]);

ax1 = axes('Position', [1 11.5 6 3]/15);
imagesc(ax1, log(edgeWeights)); 
% xticks([]); yticks([]);
title('Edge Weights');

ax2 = axes('Position', [10.5 11.5 1 3]/15);
imagesc(ax2, sum(edgeWeights,2)); 
xticks([]);
title({'Node strength', '(sum of edge weights)'});
set(ax2, 'YAxisLocation', 'right');
% colorbar;

ax3 = axes('Position', [1 9.5 6 0.5]/15);
imagesc(ax3, surfaceAreas); 
yticks([]);
title('Node Surface Area');%, 'FontWeight', 'bold');

ax4 = axes('Position', [9.5 9 3 1.5]/15);
x = surfaceAreas; y = sum(edgeWeights);
scatter(x(:), y(:), 5, 'filled'); bestfit; box on;
xticks([]); yticks([]);
xlabel('Surface Area'); ylabel('Strength');

% --------------------------------------------------------

ax5 = axes('Position', [1 5 6 3]/15);
imagesc(ax5, (surfaceAreas + surfaceAreas')); 
xticks([]); yticks([]);
title('Pairwise sum of surface area');

ax6 = axes('Position', [8 5 6 3]/15);
imagesc(ax6, (surfaceAreas .* surfaceAreas')); 
xticks([]); yticks([]);
title('Pairwise product of surface area');

ax7 = axes('Position', [1 1 6 3]/15);
x = surfaceAreas + surfaceAreas'; y = edgeWeights;
scatter(ax7, x(:), y(:), 10, 'filled'); bestfit;
xlabel('Sum of surface areas'); ylabel('Edge Weights');
xticks([]); yticks([]);
title({sprintf('Pearson Correlation = %0.2f', corr(x(:), y(:))), ...
    sprintf('Spearman Correlation = %0.2f', corr(x(:), y(:), 'Type', 'Spearman'))});

ax8 = axes('Position', [8 1 6 3]/15);
x = surfaceAreas .* surfaceAreas'; y = edgeWeights;
scatter(ax8, x(:), y(:), 10, 'filled'); bestfit;
xlabel('Product of surface areas'); ylabel('Edge Weights');
xticks([]); yticks([]);
title({sprintf('Pearson Correlation = %0.2f', corr(x(:), y(:))), ...
    sprintf('Spearman Correlation = %0.2f', corr(x(:), y(:), 'Type', 'Spearman'))});

%%
% %
% figure('Position', [400 400 600 600]);
% ax1 = axes('Position', [0.3 0.55 0.4 0.4]);
% imagesc(ax1, (surfaceAreas + surfaceAreas')); 
% xticks([]); yticks([]);
% ylabel('Pairwise sum of surface areas');
% 
% ax1 = axes('Position', [0.3 0.05 0.4 0.4]);
% imagesc(ax1, (surfaceAreas .* surfaceAreas')); 
% xticks([]); yticks([]);
% ylabel('Pairwise product of surface areas');
% 
% %
% figure('Position', [400 400 600 600]);
% tl = tiledlayout(2, 2);
% 
% nexttile();
% x = surfaceAreas; y = sum(edgeWeights);
% scatter(x(:), y(:), 'filled'); bestfit;
% xlabel('Node Surface Area'); ylabel('Node Strength'); 
% title({sprintf('Pearson Correlation = %0.4f', corr(x(:), y(:))), ...
%     sprintf('Spearman Correlation = %0.4f', corr(x(:), y(:), 'Type', 'Spearman'))});
% 
% 
% nexttile();
% 
% nexttile();
% x = surfaceAreas + surfaceAreas'; y = edgeWeights;
% scatter(x(:), y(:), 'filled');
% xlabel('Sum of surface areas');
% ylabel('Edge Weights');
% bestfit;
% title({sprintf('Pearson Correlation = %0.4f', corr(x(:), y(:))), ...
%     sprintf('Spearman Correlation = %0.4f', corr(x(:), y(:), 'Type', 'Spearman'))});
% 
% nexttile();
% x = surfaceAreas .* surfaceAreas'; y = edgeWeights;
% scatter(x(:), y(:), 'filled');
% xlabel('Product of surface areas');
% ylabel('Edge Weights');
% bestfit;
% title({sprintf('Pearson Correlation = %0.4f', corr(x(:), y(:))), ...
%     sprintf('Spearman Correlation = %0.4f', corr(x(:), y(:), 'Type', 'Spearman'))});

%% Figure: edge weight and node SA for all parcellations, pipelines, metrics, thresholds

currentDATA = ["ADJGroupWei", "ADJGroupDen", "ADJGroupConDenAltered", "ADJGroupDst"];
currentDATALabels = ["Weight", "CV", "Consistency", "Distance Bins"];
THR_LEVELS = 1:11;
parc = [1 5 3 6];



% ----

corrsAdd = zeros(length(currentDATA)*length(THR_LEVELS), length(parc)*10);
corrsMult = zeros(size(corrsAdd));

for parci = 1:length(parc)

    currentSA = mean(PARC_SA{parc(parci)}, 1);

    for metri = 1:length(currentDATA)

        DATA = eval(currentDATA(metri));

        for pipei = 1:10

            for thri = 1:length(THR_LEVELS)

                current = DATA{ORDERED_INDS{parc(parci)}(pipei), thri};
                currentX = (metri-1)*length(THR_LEVELS) + thri;
                currentY = (parci-1)*10 + pipei;

                temp = currentSA + currentSA';
                corrsAdd(currentX, currentY) = corr(temp(:), current(:));

                temp = currentSA .* currentSA';
                corrsMult(currentX, currentY) = corr(temp(:), current(:));

            end
        end
    end
end

%% corrsAdd
figure('DefaultAxesFontSize', 16, 'Position', [100 100 650 650]);
tiledlayout(4, 4, 'TileIndexing', 'columnmajor');
% cm = multi_cmap(0, 'mediumblue', 2, 'firebrick', 6);
for ii = 1:length(parc)
    for jj = 1:length(currentDATA)
        nexttile();
%         imagesc(skeww( (11*jj-10):(11*jj) , (10*ii-9):10*ii));
        imagesc(corrsAdd( (11*jj-10):(11*jj) , (10*ii-9):10*ii));
        caxis([0 1]); colormap viridis;
%         caxis([0 6]); colormap(cm);
        xticks([]);
        axis square;

        if ii == 1
            if jj == 4; ylabeldata = thr_strings_dst(1:2:end); yl = {'Number', 'of Bins'};
            else; ylabeldata = thr_strings_density(1:2:end) + "%"; yl = 'Threshold'; end

            yticks(1:2:11);
            yticklabels(ylabeldata)
            ylabel(yl);
        else
            yticks([]);
        end

    end
end

%% corrsMult
figure('DefaultAxesFontSize', 16, 'Position', [100 100 650 650]);
tiledlayout(4, 4, 'TileIndexing', 'columnmajor');
% cm = multi_cmap(0, 'mediumblue', 2, 'firebrick', 6);
for ii = 1:length(parc)
    for jj = 1:length(currentDATA)
        nexttile();
%         imagesc(skeww( (11*jj-10):(11*jj) , (10*ii-9):10*ii));
        imagesc(corrsMult( (11*jj-10):(11*jj) , (10*ii-9):10*ii));
        caxis([0 1]); colormap viridis;
%         caxis([0 6]); colormap(cm);
        xticks([]);
        axis square;

        if ii == 1
            if jj == 4; ylabeldata = thr_strings_dst(1:2:end); yl = {'Number', 'of Bins'};
            else; ylabeldata = thr_strings_density(1:2:end) + "%"; yl = 'Threshold'; end

            yticks(1:2:11);
            yticklabels(ylabeldata)
            ylabel(yl);
        else
            yticks([]);
        end

    end
end
