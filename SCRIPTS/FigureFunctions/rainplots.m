%% change density

parc = [1 5 3 6];
thr = 5;
DATA = ADJGroupConDenAltered;

lineColors = flipud([
    254	196	79;
    254	153	41;
    236	112	20;
    204	76	2;
    153	52	4;
    102	37	6]/255);

% f = @(x) sum(logical(x));
f = @(x) sum(x);

figure('DefaultAxesFontSize', 16);
tl1 = tiledlayout(1, length(parc));


t2 = [];
% for each parcellation
for kk = 1:length(parc)

    tl2 = tiledlayout(tl1, 10, 1, 'TileSpacing', 'none');
    tl2.Layout.Tile = kk;


    % get x-limits for each parcelltion
    currentDATA = DATA(ORDERED_INDS{parc(kk)}(:), :);
    currentStrengths = cellfun(f, currentDATA, 'UniformOutput', false);
    currentStrengths = [currentStrengths{:}];

    % for each pipeline
    for jj = 1:10

        nexttile(tl2);
        raincloud_plot_mean(f(DATA{ORDERED_INDS{parc(kk)}(jj), thr}), 'box_on', 1, 'color', lineColors(thr,:), 'line_width', 1);
       
        xlim([min(currentStrengths), max(currentStrengths)]);
        
        ylabel([]); yticklabels([]); yticks([]);
        yl = get(gca, 'YLim'); 
        set(gca, 'YLim', [-yl(2), yl(2)], 'XColor', [0.05 0.05 0.05]);
        
        if jj ~= 10; xlabel([]); xticklabels([]); xticks([]); end

    end


end


%% change density

parc = 5;
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
tl1 = tiledlayout(1, length(thr));


t2 = [];
% for each parcellation
for kk = 1:length(thr)

    tl2 = tiledlayout(tl1, 10, 1, 'TileSpacing', 'none');
    tl2.Layout.Tile = kk;


    % get x-limits for each parcelltion
    currentDATA = DATA(ORDERED_INDS{parc}(:), thr(kk));
    currentStrengths = cellfun(@sum, currentDATA, 'UniformOutput', false);
    currentStrengths = [currentStrengths{:}];

    % for each pipeline
    for jj = 1:10

        nexttile(tl2);
        raincloud_plot(sum(DATA{ORDERED_INDS{parc}(jj), thr(kk)}), 'box_on', 1, 'color', lineColors(kk,:), 'line_width', 1);
       
        xlim([min(currentStrengths), max(currentStrengths)]);
        
        ylabel([]); yticklabels([]); yticks([]);
        yl = get(gca, 'YLim'); 
        set(gca, 'YLim', [-yl(2), yl(2)], 'XColor', [0.05 0.05 0.05]);
        
        if jj ~= 10; xlabel([]); xticklabels([]); xticks([]); end

    end


end

