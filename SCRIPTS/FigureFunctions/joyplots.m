%% change density and parcellation

parc = [1 5 3 6];
thr = [1 3 5 7 9 11];
DATA = ADJGroupConDenAltered;

lineColors = flipud([
    254	196	79;
    254	153	41;
    236	112	20;
    204	76	2;
    153	52	4;
    102	37	6]/255);

offsets = [.00005 .00025 .00025 .00025];


figure('DefaultAxesFontSize', 16);
tl = tiledlayout(1, length(parc), 'TileSpacing', 'loose');


t2 = [];
for kk = 1:length(parc)

    % get x-limits for each parcelltion
    currentDATA = DATA(ORDERED_INDS{parc(kk)}(:), :);
    currentStrengths = cellfun(@sum, currentDATA, 'UniformOutput', false);
    currentStrengths = [currentStrengths{:}];
    e = linspace(min(currentStrengths), max(currentStrengths), 100);

    % calculate curves for each threshold (jj) and pipeline (ii)
    for jj = 1:length(thr)
        t = [];
        for ii = 1:10
            w = sum(DATA{ORDERED_INDS{parc(kk)}(ii), thr(jj)});
            t = [t; w];
        end

        for ii = 1:10
            t2(:, ii, jj) = ksdensity(t(ii,:), e)';
        end

    end

    % plot
    nexttile();
    for ii = 1:size(t2, 3)
        joyPlot(t2(:,:,ii), e, offsets(kk), 'FaceAlpha', 0, 'FaceColor', 1:10, 'LineColor', lineColors(ii,:), 'LineWidth', 2); hold on;
    end

    ylabel([]); yticks([]); yticklabels([]);
    set(gca,'view',[90 -90], 'YDir', 'reverse');

end


%% change density

parc = [5];
thr = [1 3 7 11];
DATA = ADJGroupConDenAltered;

figure;
tl = tiledlayout(1, length(thr), 'TileSpacing', 'loose');
t2 = [];

for jj = 1:length(thr)
    nexttile();
    t = [];
    for ii = 1:10
        w = sum(DATA{ORDERED_INDS{parc}(ii), thr(jj)});
        t = [t; w];
    end

    [~,e] = histcounts(t(:), 100);


    for ii = 1:10
        t2(:, ii, jj) = ksdensity(t(ii,:), e)';
    end

    joyPlot(atan(t2(:,:,jj)), e, .00010, 'FaceColor', 1:10);
    colormap plasma;

end

%%
figure;
lineColors = flipud([
    254	196	79;
    254	153	41;
    236	112	20;
    204	76	2;
    153	52	4;
    102	37	6]/255);

for ii = 1:size(t2, 3)
    joyPlot(t2(:,:,ii), e, .00010, 'FaceAlpha', 0, 'FaceColor', 1:10, 'LineColor', lineColors(ii,:), 'LineWidth', 2); hold on;
end


%% change density and parcellation

parc = [1 5 3 6];
thr = [1 3 5 7 9 11];
DATA = ADJGroupConDenAltered;

lineColors = flipud([
    254	196	79;
    254	153	41;
    236	112	20;
    204	76	2;
    153	52	4;
    102	37	6]/255);

offsets = [.00005 .00025 .00025 .00025];


figure;
tl = tiledlayout(1, length(parc), 'TileSpacing', 'loose');


t2 = [];
for kk = 1:length(parc)

    % get x-limits for each parcelltion
    currentDATA = DATA(ORDERED_INDS{parc(kk)}(:), :);
    currentStrengths = cellfun(@sum, currentDATA, 'UniformOutput', false);
    currentStrengths = [currentStrengths{:}];
    e = linspace(min(currentStrengths), max(currentStrengths), 100);

    % calculate curves for each threshold (jj) and pipeline (ii)
    for jj = 1:length(thr)
        t = [];
        for ii = 1:10
            w = sum(DATA{ORDERED_INDS{parc(kk)}(ii), thr(jj)});
            t = [t; w];
        end

        for ii = 1:10
            t2(:, ii, jj) = ksdensity(t(ii,:), e)';
        end

    end

    % plot
    nexttile();
    for ii = 1:size(t2, 3)
        joyPlot(t2(:,:,ii), e, offsets(kk), 'FaceAlpha', 0, 'FaceColor', 1:10, 'LineColor', lineColors(ii,:), 'LineWidth', 2); hold on;
    end



end


%% change density
%
% parc = [5];
% thr = [1 3 7 11];
% DATA = ADJGroupConDenAltered;
%
% figure;
% tl = tiledlayout(1, length(thr), 'TileSpacing', 'loose');
% t2 = [];
%
% for jj = 1:length(thr)
%     nexttile();
%     t = [];
%     for ii = 1:10
%         w = sum(DATA{ORDERED_INDS{parc}(ii), thr(jj)});
%         t = [t; w];
%     end
%
%     [~,e] = histcounts(t(:), 100);
%
%
%     for ii = 1:10
%         t2(:, ii, jj) = ksdensity(t(ii,:), e)';
%     end
%
%     joyPlot(t2(:,:,jj), e, .00010, 'FaceColor', 1:10);
%     colormap plasma;
%
% end
%
% %%
% figure;
% lineColors = flipud([
%     254	196	79;
%     254	153	41;
%     236	112	20;
%     204	76	2;
%     153	52	4;
%     102	37	6]/255);
%
% for ii = 1:size(t2, 3)
%     joyPlot(t2(:,:,ii), e, .00010, 'FaceAlpha', 0, 'FaceColor', 1:10, 'LineColor', lineColors(ii,:), 'LineWidth', 2); hold on;
% end
%



%% change density and parcellation

parc = [1 5 3 6];
thr = [7];
DATA = ADJGroupConDenAltered;

lineColors = flipud([
    254	196	79;
    254	153	41;
    236	112	20;
    204	76	2;
    153	52	4;
    102	37	6]/255);

offsets = [.00001 .00005 .00005 .00005];


figure;
tl = tiledlayout(1, length(parc), 'TileSpacing', 'loose');


t2 = [];
for kk = 1:length(parc)

    % get x-limits for each parcelltion
    currentDATA = DATA(ORDERED_INDS{parc(kk)}(:), :);
    currentStrengths = cellfun(@sum, currentDATA, 'UniformOutput', false);
    currentStrengths = [currentStrengths{:}];
    e = linspace(min(currentStrengths), max(currentStrengths), 100);

    % calculate curves for each threshold (jj) and pipeline (ii)
    for jj = 1:length(thr)
        t = [];
        for ii = 1:10
            w = sum(DATA{ORDERED_INDS{parc(kk)}(ii), thr(jj)});
            t = [t; w];
        end

        for ii = 1:10
            t2(:, ii, jj) = ksdensity(t(ii,:), e)';
        end

    end

    % plot
    nexttile();
    for ii = 1:size(t2, 3)
        joyPlot(t2(:,:,ii), e, offsets(kk), 'FaceAlpha', 0, 'FaceColor', 1:10, 'LineColor', lineColors(ii,:), 'LineWidth', 2); hold on;
    end



end


%% change density
%
% parc = [5];
% thr = [1 3 7 11];
% DATA = ADJGroupConDenAltered;
%
% figure;
% tl = tiledlayout(1, length(thr), 'TileSpacing', 'loose');
% t2 = [];
%
% for jj = 1:length(thr)
%     nexttile();
%     t = [];
%     for ii = 1:10
%         w = sum(DATA{ORDERED_INDS{parc}(ii), thr(jj)});
%         t = [t; w];
%     end
%
%     [~,e] = histcounts(t(:), 100);
%
%
%     for ii = 1:10
%         t2(:, ii, jj) = ksdensity(t(ii,:), e)';
%     end
%
%     joyPlot(t2(:,:,jj), e, .00010, 'FaceColor', 1:10);
%     colormap plasma;
%
% end
%
% %%
% figure;
% lineColors = flipud([
%     254	196	79;
%     254	153	41;
%     236	112	20;
%     204	76	2;
%     153	52	4;
%     102	37	6]/255);
%
% for ii = 1:size(t2, 3)
%     joyPlot(t2(:,:,ii), e, .00010, 'FaceAlpha', 0, 'FaceColor', 1:10, 'LineColor', lineColors(ii,:), 'LineWidth', 2); hold on;
% end
%





