%% Figure S3a: skewness and kurtosis for strength distributions for all parcellations, pipelines, metrics, thresholds

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

corrsw = zeros(size(skeww));
corrsd = zeros(size(skeww));

skewdr = zeros(size(skeww));
kurtdr = zeros(size(skeww));

residualCorr = zeros(size(skeww));
testResiduals = zeros(size(skeww));

perw = zeros(size(skeww)); % Jordanova Petkova 2017
perd = zeros(size(skeww));

for parci = 1:length(parc)

    currentSA = mean(PARC_SA{parc(parci)}, 1);

    for metri = 1:length(currentDATA)

        DATA = eval(currentDATA(metri));

        for pipei = 1:10

            for thri = 1:length(THR_LEVELS)

                current = DATA{ORDERED_INDS{parc(parci)}(pipei), thri};
                w = sum(current);
                d = sum(logical(current));

                %                 p = polyfit(currentSA, w,1);
                %                 wHat = polyval(p,currentSA);
                %                 wr = w-wHat; % weight_residuals
                wr = calcResiduals(currentSA, w);
                dr = calcResiduals(currentSA, d);

                currentX = (metri-1)*length(THR_LEVELS) + thri;
                currentY = (parci-1)*10 + pipei;

                skeww(currentX, currentY) = skewness(w);
                skewwr(currentX, currentY) = skewness(wr);

                skewd(currentX, currentY) = skewness(d);
                skewdr(currentX, currentY) = skewness(dr);

                kurtw(currentX, currentY) = kurtosis(w) - 3;
                kurtwr(currentX, currentY) = kurtosis(wr) - 3;

                kurtd(currentX, currentY) = kurtosis(d) - 3;
                kurtdr(currentX, currentY) = kurtosis(dr) - 3;

                corrsw(currentX, currentY) = corr(currentSA', w', 'Type', 'Pearson');
                corrsd(currentX, currentY) = corr(currentSA', d', 'Type', 'Pearson');
                residualCorr(currentX, currentY) = corr(w', wr');

                testResiduals(currentX, currentY) = corr(wr', currentSA');

                perw(currentX, currentY) = pExtremeRight(w);
                pmrw(currentX, currentY) = pMildRight(w);
                perd(currentX, currentY) = pExtremeRight(d);
                pmrd(currentX, currentY) = pMildRight(d);
                perwr(currentX, currentY) = pExtremeRight(wr);
                perdr(currentX, currentY) = pExtremeRight(dr);

            end
        end
    end
end
%% TailIndices
% figure;
% subplot(2, 2, 1); imagesc(pmrw);
% caxis([0, 0.0389]); cm = multi_cmap(0, 'mediumblue', 0.0389, 'firebrick', 0.0389*3); colormap(cm);
% title('Weighted, 1.5*IQR');
% 
% subplot(2, 2, 3); imagesc(pmrd);
% caxis([0, 0.0389]); cm = multi_cmap(0, 'mediumblue', 0.0389, 'firebrick', 0.0389*3); colormap(cm);
% title('Binarised, 1.5*IQR');
% 
% subplot(2, 2, 2); imagesc(perw);
% caxis([0, 3/108]); cm = multi_cmap(0, 'mediumblue', 1/108, 'firebrick', 3/108); colormap(cm);
% title('Weighted, 3*IQR');
% 
% subplot(2, 2, 4); imagesc(perd);
% caxis([0, 3/108]); cm = multi_cmap(0, 'mediumblue', 1/108, 'firebrick', 3/108); colormap(cm);
% title('Binarised, 3*IQR');
figure('DefaultAxesFontSize', 16, 'Position', [100 100 650 650]);
tiledlayout(4, 4, 'TileIndexing', 'columnmajor');
cm = multi_cmap(0, 'mediumblue', 1/108, 'firebrick', 3/108);
for ii = 1:length(parc)
    for jj = 1:length(currentDATA)
        nexttile();
        %         imagesc(skeww( (11*jj-10):(11*jj) , (10*ii-9):10*ii));
        % -----------------------------------------------------------------
        imagesc(perw( (11*jj-10):(11*jj) , (10*ii-9):10*ii));
        % -----------------------------------------------------------------
        %         caxis([0 2]); colormap viridis;
        caxis([0 3/108]); colormap(cm);
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


%% skewness
figure('DefaultAxesFontSize', 16, 'Position', [100 100 650 650]);
tiledlayout(4, 4, 'TileIndexing', 'columnmajor');
cm = multi_cmap(0, 'mediumblue', 2, 'firebrick', 6);
for ii = 1:length(parc)
    for jj = 1:length(currentDATA)
        nexttile();
        %         imagesc(skeww( (11*jj-10):(11*jj) , (10*ii-9):10*ii));
        imagesc(skewd( (11*jj-10):(11*jj) , (10*ii-9):10*ii));
        %         caxis([0 2]); colormap viridis;
        caxis([0 6]); colormap(cm);
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

%% kurtosis

figure('DefaultAxesFontSize', 16, 'Position', [100 100 650 650]);
tiledlayout(4, 4, 'TileIndexing', 'columnmajor');
cm = multi_cmap(0, 'mediumblue', 6, 'firebrick', 18);
for ii = 1:length(parc)
    for jj = 1:length(currentDATA)
        nexttile();
        %         imagesc(kurtw( (11*jj-10):(11*jj) , (10*ii-9):10*ii));
        imagesc(kurtw( (11*jj-10):(11*jj) , (10*ii-9):10*ii));
        %         caxis([0 6]); colormap viridis;
        caxis([0 18]);  colormap(cm);
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

%% correlation

figure('DefaultAxesFontSize', 16, 'Position', [100 100 650 650]);
tiledlayout(4, 4, 'TileIndexing', 'columnmajor');

for ii = 1:length(parc)
    for jj = 1:length(currentDATA)
        nexttile();
        %         imagesc(corrsw( (11*jj-10):(11*jj) , (10*ii-9):10*ii));
        imagesc(corrsd( (11*jj-10):(11*jj) , (10*ii-9):10*ii));
        caxis([0 1]); colormap viridis;
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

%% correlation with residuals

figure('DefaultAxesFontSize', 16, 'Position', [100 100 650 650]);
tiledlayout(4, 4, 'TileIndexing', 'columnmajor');

for ii = 1:length(parc)
    for jj = 1:length(currentDATA)
        nexttile();
        imagesc(residualCorr( (11*jj-10):(11*jj) , (10*ii-9):10*ii));
        caxis([0 1]); colormap viridis;
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

%% skewness and kurtosis compared to original

defColors = [0, 0.6470, 0.9410;
    0.8500, 0.5250, 0.2980;
    0.9290, 0.3940, 0.3250;
    0.6940, 0.3840, 0.7560;
    0.4660, 0.6740, 0.1880;
    0.3010, 0.7450, 0.9330;
    0.6350, 0.0780, 0.1840];

marker = {'o', 's', 'd', '^', 'p'};
markerEdgeColor = {[0.5 0.5 0.5], 'k'};

figure('DefaultAxesFontSize', 14, 'Position', [100 100 600 750]); hold on;

tl = tiledlayout(1, 2);

for kk = 1:2
    tl2 = tiledlayout(tl, 3, 1, 'TileSpacing', 'loose');
    tl2.Layout.Tile = kk;

    nexttile(tl2); % skewness
    if kk == 1; dataA = skeww; dataB = skewwr;
    elseif kk == 2; dataA = skewd; dataB = skewdr; end
    for ii = 1:10; for jj = 1:4; hold on;
            s = scatter( dataA( :,10*(jj-1)+ii ) , dataB( :,10*(jj-1)+ii ), ...
                50, defColors(jj,:), 'filled', 'MarkerEdgeColor','k');%, ...
            %             marker{mod(ii-1,5)+1}, 'MarkerEdgeColor', markerEdgeColor{(ii > 5.5)+1}, 'LineWidth', 1);
            s.MarkerFaceAlpha = 0.35;
    end; end
axis square; box on; same_axes;
xlabel('Original Skewness');
ylabel({'Skewness',' of Residuals'});
plot(xlim, ylim, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');


    nexttile(tl2); % righttailedness
    if kk == 1; dataA = perw; dataB = perwr;
    elseif kk == 2; dataA = perd; dataB = perdr; end
    for ii = 1:10; for jj = 1:4; hold on;
            s = scatter( dataA( :,10*(jj-1)+ii ) , dataB( :,10*(jj-1)+ii ), ...
                50, defColors(jj,:), 'filled', 'MarkerEdgeColor','k');%, ...
            %             marker{mod(ii-1,5)+1}, 'MarkerEdgeColor', markerEdgeColor{(ii > 5.5)+1}, 'LineWidth', 1);
            s.MarkerFaceAlpha = 0.35;
    end; end
axis square; box on; same_axes;
xlabel('Original Right-tailedness');
ylabel({'Right-tailedness',' of Residuals'});
plot(xlim, ylim, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');


    nexttile(tl2); % kurtosis
    if kk == 1; dataC = kurtw; dataD = kurtwr;
    elseif kk == 2; dataC = kurtd; dataD = kurtdr; end
    for ii = 1:10; for jj = 1:4; hold on;
            s = scatter( dataC( :,10*(jj-1)+ii ) , dataD( :,10*(jj-1)+ii ), ...
                50, defColors(jj,:), 'filled', 'MarkerEdgeColor','k');%, ...
        %             marker{mod(ii-1,5)+1}, 'MarkerEdgeColor', markerEdgeColor{(ii > 5.5)+1}, 'LineWidth', 1);
            s.MarkerFaceAlpha = 0.35;
    end; end
axis square; box on; same_axes;
xlabel('Original Kurtosis');
ylabel({'Kurtosis',' of Residuals'});
plot(xlim, ylim, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');


if kk == 1; title(tl2, 'Weighted', 'FontSize', 20, 'FontWeight', 'bold');
elseif kk == 2; title(tl2, 'Unweighted', 'FontSize', 20, 'FontWeight', 'bold'); end

end

% l = legend(parc_name2{[1 5 3 6]}, 'Location', 'eastoutside');
legend(parc_name2{[1 5 3 6]}, 'Location', 'eastoutside');
% l.Title.String = "Parcellation";
% l.Title.Visible = 'on';

%%

for ii = 1:10
    for jj = 1:4

        s = scatter( skeww( :,10*(jj-1)+ii ) , skewwr( :,10*(jj-1)+ii ), ...
            100, defColors(jj,:), 'filled', 'MarkerEdgeColor','k');%, ...
        %             marker{mod(ii-1,5)+1}, 'MarkerEdgeColor', markerEdgeColor{(ii > 5.5)+1}, 'LineWidth', 1);
        s.MarkerFaceAlpha = 0.5;

    end
end

axis square; box on; same_axes;
xlabel('Original Skewness');
ylabel('Skewness of Residuals');

plot(xlim, ylim, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');

l = legend(parc_name2{[1 5 3 6]}, 'Location', 'eastoutside');
l.Title.String = "Parcellation";

% ---

figure('DefaultAxesFontSize', 20, 'Position', [100 100 640 400]); hold on;

for ii = 1:10
    for jj = 1:4

        s = scatter( kurtw( :,10*(jj-1)+ii ) , kurtwr( :,10*(jj-1)+ii ), ...
            100, defColors(jj,:), 'filled', 'MarkerEdgeColor','k');%, ...
        %             marker{mod(ii-1,5)+1}, 'MarkerEdgeColor', markerEdgeColor{(ii > 5.5)+1}, 'LineWidth', 1);
        s.MarkerFaceAlpha = 0.5;

    end
end

axis square; box on; same_axes;
xlabel('Original Kurtosis');
ylabel('Kurtosis of Residuals');

plot(xlim, ylim, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':');

l = legend(parc_name2{[1 5 3 6]}, 'Location', 'eastoutside');
l.Title.String = "Parcellation";

%%
%
%
% %%
%
% figure('Position', [-1.3206e+03 221.8000 1.1808e+03 544], 'DefaultAxesFontSize', 18);
% tiledlayout(2, 2, 'TileSpacing','compact', 'Padding', 'compact')
%
% nexttile();
% imagesc(skeww);
% title({'Skewness', 'Strength'});
% yticks([6, 17, 28, 39, 50]); caxis([0 2]);
% yticklabels(currentDATALabels);
%
% nexttile();
% imagesc(skewwr);
% title({'Skewness', 'Strength Residuals'});
% yticks([6, 17, 28, 39, 50]); caxis([0 2]);
% yticklabels([]); colorbar;
%
% nexttile();
% imagesc(kurtw-3);
% title({'Excess Kurtosis', 'Strength'});
% yticks([6, 17, 28, 39, 50]); caxis([0 6]);
% yticklabels(currentDATALabels);
%
% nexttile();
% imagesc(kurtwr-3);
% title({'Excess Kurtosis', 'Strength Residuals'});
% yticks([6, 17, 28, 39, 50]); caxis([0 6]);
% yticklabels([]); colorbar;
%
% for ii = 1:4
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
%     %     axis square;  %caxis([0 1]);
%     xticks([5.5, 15.5, 25.5, 35.5]);
%     xticklabels(["82", "S220", "HCP", "S520"]);
%
% end
%
% scfh(1000); scfw(1000);
% % print('v3_4a', '-dpng', '-r600');
%
%
%
% %%
% currentDATA = ["ADJGroupWei", "ADJGroupDen", "ADJGroupConDenAltered", "ADJGroupDst"];
% currentDATALabels = ["Weight", "CV", "Consistency", "Distance Bins"];
% THR_LEVELS = 1:11;
% parc = [1 5 3 6];
%
%
%
% for parci = 1:length(parc)
%
%     currentSA = mean(PARC_SA{parc(parci)}, 1);
%
%     for metri = 1:length(currentDATA)
%
%         DATA = eval(currentDATA(metri));
%
%         for pipei = 1:10
%
%             for thri = 1:length(THR_LEVELS)
%
%
%                 current = DATA{ORDERED_INDS{parc(parci)}(pipei), thri};
%                 w = sum(current);
%
%
%                 [~,~,timr(... % tail index mild right tail
%                     (metri-1)*length(THR_LEVELS) + thri, ...
%                     (parci-1)*10 + pipei                    ),...
%                     tier(... % tail index extreme right tail
%                     (metri-1)*length(THR_LEVELS) + thri, ...
%                     (parci-1)*10 + pipei                    )] = getFenceIndices(w);
%
%             end
%         end
%     end
% end
%
%
% %%
%
%
% figure('Position', [-1.3206e+03 221.8000 1.1808e+03 544], 'DefaultAxesFontSize', 18);
% tiledlayout(1, 2, 'TileSpacing','compact', 'Padding', 'compact')
%
% nexttile();
% imagesc(timr);
% title({'Intermediate Right Tail', 'Strength'});
% yticks([6, 17, 28, 39, 50]);
% yticklabels(currentDATALabels); colorbar; caxis([0 0.0339]);
%
% nexttile();
% imagesc(tier);
% title({'Extreme Right Tail', 'Strength'});
% yticks([6, 17, 28, 39, 50]); yticklabels([]);colorbar;
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
%     axis square;
%     xticks([5.5, 15.5, 25.5, 35.5]);
%     xticklabels(["82", "S220", "HCP", "S520"]);
%
% end
%
%
% scfh(500); scfw(1000);
% print('v3_4b', '-dpng', '-r600');
%
% %% 4c example of kurtosis and skewness
%
% parc = 3;
% pipeline = 6;
% currentDATA = "ADJGroupConDen";
%
% figure;
% tiledlayout(1, 11);
% DATA = eval(currentDATA);
%
% clear ax;
% for ii = 1:11
%
%     ax(ii) = nexttile();
%
%
%     current = DATA{ORDERED_INDS{parc}(pipeline), THR_LEVELS(ii)};
%     [epdf, xi] = ksdensity(sum(current~=0), 'Function', 'pdf');
%     plot(xi, epdf, 'LineWidth', 2);
%
%
% end
%
% linkaxes(ax);