%% Figure 2a: plots example of changes along each of the axes
% with defaults specified by the user

% user to adjust defaults
defParc = 2;
defPipe = 5;
defGroup = 2;
defLevel = 4;

% ---
parcs = [1 5 3 6];
parcLabels = parc_name2;

pipelines = 1:10;
pipeLabels = pipeline_titles;

groupNames = {ADJGroupWei, ADJGroupDen, ADJGroupConDenAltered, ADJGroupDst};
groupLabels = {'Weight', 'CV', 'Consistency', 'Distance Dependence'};

thrLevels = 1:2:11;
thrLabels = cellfun(@(x) x + "%", thr_strings_density, 'UniformOutput', false);

verts = [lh_inflated_verts; rh_inflated_verts];
nodeLocs = {[lh_aparc; rh_aparc], [lh_rand200; rh_rand200], [lh_HCPMMP1; rh_HCPMMP1],  [lh_rand500; rh_rand500], [Scha_parcs.lh_scha200; Scha_parcs.rh_scha200], [Scha_parcs.lh_scha500; Scha_parcs.rh_scha500]};


%% ---

figure('Position', [1200 200 820 840]);
t1 = tiledlayout(5, 1, 'TileSpacing', 'none', 'Padding', 'none');
tiles = {1, [2, 3], 4, 5};

for ii = 1:4

    % get n and also set the parameters to be used for plotting
    switch ii
        case 1 % Parcellation
            n = length(parcs); spec = defParc;
            yt = 'Parcellation'; pt = parcLabels(parcs);
        case 2 % Pipeline
            n = length(pipelines); spec = defPipe;
            yt = 'Pipeline'; pt = pipeLabels;
        case 3 % Group
            n = length(groupNames); spec = defGroup;
            yt = {'Group', 'Reconstruction'}; pt = groupLabels;
        case 4 % Level
            n = length(thrLevels); spec = defLevel;
            yt = 'Threshold'; pt = thrLabels(thrLevels);
    end

    currParc = repmat(parcs(defParc), 1, n);
    currPipe = repmat(defPipe, 1, n);
    currGroup = repmat(defGroup, 1, n);
    currLevel = repmat(thrLevels(defLevel), 1, n);

    switch ii
        case 1 % Parcellation
            currParc = parcs;
        case 2 % Pipeline
            currPipe = pipelines;
        case 3 % Group
            currGroup = 1:n;
        case 4 % Level
            currLevel = thrLevels;
    end


    % do plotting (needs locations, SA, and edge weights)
    t2 = tiledlayout(t1, tern(ii==2, 2, 1), n/tern(ii==2, 2, 1), 'TileSpacing', 'none', 'Padding', 'none');
    t2.Layout.Tile = tiles{ii};
    t2.Layout.TileSpan = [length(tiles{ii}), 1];

    clear ax;
    for jj = 1:n
        ax(jj) = nexttile(t2, jj);

        % node locations
        rois = joinParcs(nodeLocs{currParc(jj)});
%         nodeLocations = zeros(max(rois), 3);
%         for temp = 1:max(rois)
%             nodeLocations(temp,:) = mean(verts(rois == temp, :), 1);
%         end

        % node SA
        currentSA = mean(PARC_SA{currParc(jj)}, 1);

        % node edge weights
        A = groupNames{currGroup(jj)};
        A = A{ORDERED_INDS{currParc(jj)}(currPipe(jj)), currLevel(jj)};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%         w = sum(A); 
        w = sum(logical(A));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        [~,worder] = sort(w, 'ascend');

        if jj == 1
            ylabel(yt, 'FontSize', 16, 'FontWeight', 'bold');
            ax.YLabel.Visible = 'on';
        end

        t3 = tiledlayout(t2, 2, 1, 'TileSpacing', 'none', 'Padding', 'none');
        t3.Layout.Tile = jj;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        % for titles:
% title(split(pt{jj}, ["/", " "]));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% for images: 
        nexttile(t3);
        brainplot(verts(1:end/2,:), faces, rois(1:end/2,:), w(1:end/2), 0, plasma);
%         ylabel(split(pt{jj}, '/')); set(get(gca, 'YLabel'), 'Rotation', 0);
        colorbar off; 

        nexttile(t3);
        brainplot(verts(1:end/2,:), faces, rois(1:end/2,:), w(1:end/2), 0, plasma);
        view([90 0]); colorbar off;
        
       
        axis(ax(jj), 'off');
%         scatter(nodeLocations(worder,1), nodeLocations(worder,2), atan(currentSA(worder)/800)*50, ...
%             (w(worder)), 'filled', 'MarkerEdgeColor', 'black');

        axis equal; xlim([-75 75]); ylim([-75 100]); axis off;

%         if jj == spec
%             hold on; rectangle('Position', [-75 -75 150 165]); hold off;
%         end
%         if jj == 1
%             ylabel(yt, 'FontSize', 16, 'FontWeight', 'bold');
%             ax.YLabel.Visible = 'on';
%         end
%         title(pt{jj});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    end

end


% %% Figure 2b: examples of degree distributions (all combinations)
% parcs = [1 5 3];
% parcLabels = parc_name;
% 
% pipelines = 1:10;
% pipeLabels = pipeline_titles;
% 
% groupNames = {ADJGroupWei, ADJGroupDen, ADJGroupConDenAltered, ADJGroupDst};
% groupLabels = {'Weight', 'CV', 'Consistency', 'Distance Dependence'};
% 
% thrLevels = [1 3 5 7 11];
% thrLabels = cellfun(@(x) x + "%", thr_strings_density, 'UniformOutput', false);
% 
% colors = turbo(length(thrLevels));
% 
% 
% %%---
% 
% figure("Position", [-1500 200 1427 432]);
% t1 = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
% 
% 
% for ii = 1:length(groupNames)
% 
%     t2 = tiledlayout(t1, length(parcs), length(pipelines), 'TileSpacing', 'none', 'Padding', 'none');
%     t2.Layout.Tile = ii;
% 
%     currentGroup = groupNames{ii};
% 
%     for jj = 1:length(parcs)
% 
%         for kk = 1:length(pipelines)
% 
%             nexttile(t2);
% 
%             for ll = 1:length(thrLevels)
% 
%                 current = currentGroup{ORDERED_INDS{parcs(jj)}(pipelines(kk)), thrLevels(ll)};
%                 d = sum(current);
% 
%                 [epdf, xi] = ksdensity(d, 'Function', 'pdf');
%                 plot(xi, epdf, 'LineWidth', 2, 'LineStyle','-', 'Color', colors(ll,:));
%                 hold on; axis square; xticks([]); yticks([]);
% 
%             end
% 
%         end
%     end
% 
% end





%% Figure 2b_v2 (in the style of 2a)
% with defaults specified by the user
thrLevels = [1 3 5 7 9 11];

% user to adjust defaults
defParc = 2;
defPipe = 5;
defGroup = 2;

colors = flipud([
254	196	79;
254	153	41;
236	112	20;
204	76	2;
153	52	4;
102	37	6]/255);
% colors = colors(end:-1:1,:);

% ---
parcs = [1 5 3 6];
parcLabels = parc_name2;

pipelines = 1:10;
pipeLabels = pipeline_titles;

groupNames = {ADJGroupWei, ADJGroupDen, ADJGroupConDenAltered, ADJGroupDst};
groupLabels = {'Weight', 'CV', 'Consistency', 'Distance Dependence'};

% thrLevels = [1 3 5 7 9 11];
thrLabels = cellfun(@(x) x + "%", thr_strings_density, 'UniformOutput', false);

% verts = [lh_verts; rh_verts];
% nodeLocs = {[lh_aparc; rh_aparc], [lh_rand200; rh_rand200], [lh_HCPMMP1; rh_HCPMMP1],  [lh_rand500; rh_rand500], [Scha_parcs.lh_scha200; Scha_parcs.rh_scha200], [Scha_parcs.lh_scha500; Scha_parcs.rh_scha500]};


% %% ---

figure('Position', [-1351 200 1280 480]);
t1 = tiledlayout(4, 1, 'TileSpacing', 'loose', 'Padding', 'none');
tiles = {1, [2, 3], 4};

for ii = 1:3

    % get n and also set the parameters to be used for plotting
    switch ii
        case 1 % Parcellation
            n = length(parcs); spec = defParc;
            yt = 'Parcellation'; pt = parcLabels(parcs);
        case 2 % Pipeline
            n = length(pipelines); spec = defPipe;
            yt = 'Pipeline'; pt = pipeLabels;
        case 3 % Group
            n = length(groupNames); spec = defGroup;
            yt = {'Group', 'Reconstruction'}; pt = groupLabels;
    end

    currParc = repmat(parcs(defParc), 1, n);
    currPipe = repmat(defPipe, 1, n);
    currGroup = repmat(defGroup, 1, n);

    switch ii
        case 1 % Parcellation
            currParc = parcs;
        case 2 % Pipeline
            currPipe = pipelines;
        case 3 % Group
            currGroup = 1:n;
    end
    

    % do plotting (needs locations, SA, and edge weights)
    t2 = tiledlayout(t1, tern(ii==2, 2, 1), n/tern(ii==2, 2, 1), 'TileSpacing', 'compact', 'Padding', 'none');
    t2.Layout.Tile = tiles{ii};
    t2.Layout.TileSpan = [length(tiles{ii}), 1];

    clear ax;
    for jj = 1:n
        ax(jj) = nexttile(t2, jj);

        % node locations
%         rois = joinParcs(nodeLocs{currParc(jj)});
%         nodeLocations = zeros(max(rois), 3);
%         for temp = 1:max(rois)
%             nodeLocations(temp,:) = mean(verts(rois == temp, :), 1);
%         end

        % node SA
        currentSA = mean(PARC_SA{currParc(jj)}, 1);

        % node edge weights
        currentGroup = groupNames{currGroup(jj)};

        for kk = 1:length(thrLevels)
            A = currentGroup{ORDERED_INDS{currParc(jj)}(currPipe(jj)), thrLevels(kk)};
            d = sum(A);

            [epdf, xi] = ksdensity(d, 'Function', 'pdf');
            plot(xi, epdf, 'LineWidth', 3, 'LineStyle','-', 'Color', colors(kk,:));
            hold on; axis square; xticks([]); yticks([]);
        end

%         if jj == spec
%             ax(jj).LineWidth = 3;
%         end
%         if jj == 1
%             ylabel(yt);
%             ax.YLabel.Visible = 'on';
%         end
        ylabel(split(pt{jj}, ["/", " "]), 'FontSize', 12, 'FontWeight', 'bold'); 
        set(get(gca, 'YLabel'), 'Rotation', 0, 'VerticalAlignment','middle', 'HorizontalAlignment','right');

    end

end

%% Make Legend
figure; hold on; 
for ii = 1:length(thrLevels); plot([0 0], 'color', colors(ii,:), 'LineWidth', 5); end; 
l = legend(thrLabels(thrLevels), 'Location', 'eastoutside', 'FontSize', 30); 
l.Title.String = 'Threshold'; 
axis off
set(gca, 'Position', [0 0 0 1]);
scfh(375); scfw(212);
% set(l, 'Position', [0 0 1 1]);


