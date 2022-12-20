%% change density

parc = [5];
thr = [1 3 5 7 9 11];
DATA = ADJGroupConDenAltered;

lineColors = flipud([
    254	196	79;
    254	153	41;
    236	112	20;
    204	76	2;
    153	52	4;
    102	37	6]/255);

figure('DefaultAxesFontSize', 16);
tl = tiledlayout(10, 3, 'TileSpacing', 'none');

% get x-limits for the parcelltion
currentDATA = DATA(ORDERED_INDS{parc}(:), :);
currentStrengths = cellfun(@sum, currentDATA, 'UniformOutput', false);
currentStrengths = [currentStrengths{:}];
e = linspace(min(currentStrengths), max(currentStrengths), 100);


for ii = 1:10

    nexttile();


    t2 = [];
    % calculate curves for each threshold (jj)
    for jj = 1:length(thr)

        w = sum(DATA{ORDERED_INDS{parc}(ii), thr(jj)});


        t2(:, jj) = ksdensity(w, e)';

        hold on;
        plot(e, t2);



    end

    set(gca, 'View', [-90 90]);


    nexttile([1 2]);

    t2 = [];
    for jj = 1:length(thr)
        t2(:,jj) = sum(DATA{ORDERED_INDS{parc}(ii), thr(jj)})';

    end

    boxplot(t2);

end


%% change density

parc = [6];
thr = [1 3 5 7 9 11];
DATA = ADJGroupConDenAltered;

lineColors = flipud([
    254	196	79;
    254	153	41;
    236	112	20;
    204	76	2;
    153	52	4;
    102	37	6]/255);

SCALING = 1e4;

figure;
tl = tiledlayout(3, 10, 'TileSpacing', 'none', 'TileIndexing', 'columnmajor');

% get x-limits for the parcelltion
currentDATA = DATA(ORDERED_INDS{parc}(:), :);
currentStrengths = cellfun(@sum, currentDATA, 'UniformOutput', false);
currentStrengths = [currentStrengths{:}]/SCALING;
e = linspace(min(currentStrengths), max(currentStrengths), 100);


for ii = 1:10

    nexttile();
    t2 = [];
    % calculate curves for each threshold (jj)
    for jj = 1:length(thr)
        w = sum(DATA{ORDERED_INDS{parc}(ii), thr(jj)})/SCALING;
        [t2, e] = ksdensity(w);
        hold on;
        plot(e, t2, 'LineWidth', 2, 'Color', lineColors(jj,:));
    end
%     xlim([min(e), max(e)]);
    xlabel([]); xticks([]);
    ylabel([]); yticks([]);
    box on;


    nexttile([2 1]);
    t2 = [];
    for jj = 1:length(thr)
        t2(:,jj) = sum(DATA{ORDERED_INDS{parc}(ii), thr(jj)})'/SCALING;
    end
    bp = boxplot(fliplr(t2), 'Orientation', 'horizontal', 'Color', flipud(lineColors), 'Symbol', '|');
    for ibp = 1:7; set(bp(ibp,:), 'LineWidth', 1); end
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');

    if ii == 1
        yticks(1:length(thr));
        yticklabels(thr_strings_density(fliplr(thr)) + "%");
        ylabel('Density', 'FontSize', 18);
    else
        yticks([]);
    end
    

end

xlabel(tl, 'Node Strength (\times10^4)', 'FontSize', 18);
title(tl, parc_name2(parc));
scfw(800); scfh(300);