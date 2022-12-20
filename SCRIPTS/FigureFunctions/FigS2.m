%% Fig S2a: Hub locations as a function of parcellation and pipeline

verts = lh_inflated_verts;
faces = lh_faces;
currentDATA = 'ADJGroupDen';
THR_LEVELS = 5;
parc = [1 5 3 6];
roisAll = { [lh_aparc], ...
    [Scha_parcs.lh_scha200], ...
    [lh_HCPMMP1], [Scha_parcs.lh_scha500]};

% ----

clear ax;
ax = gobjects(10, length(parc));

DATA = eval(currentDATA);

f = figure('Position', [500, 450, 160*length(parc) , 1100], 'DefaultAxesFontSize', 30);
t1 = tiledlayout(1, length(parc), "TileSpacing", 'none', 'Padding', 'none');



for jj = 1:length(parc)



    % re-index rois in right hemisphere
    rois = roisAll{jj};

    nodeLocations = zeros(max(rois), 3);
    for ii = 1:max(rois)
        nodeLocations(ii,:) = mean(verts(rois == ii, :), 1);
    end

    %     currentSA = mean(PARC_SA{parc(jj)}, 1);


    t2 = tiledlayout(t1, 10, 1, 'TileSpacing', 'none', 'Padding', 'compact');
    t2.Layout.Tile = jj;

    for ii = 1:10


        current = DATA{ORDERED_INDS{parc(jj)}(ii), THR_LEVELS};
        w = sum(current);
        [~, worder] = sort(w, 'ascend');

        ax(ii,jj) = nexttile(t2); axis(ax(ii,jj), 'off');

        t3 = tiledlayout(t2, 1, 2, 'TileSpacing', 'none', 'Padding', 'compact');
        t3.Layout.Tile = ii;


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

        %         if ii == 1
        %             ylabel(ax(ii,jj), parc_name{parc(jj)});
        %             ax(ii,jj).YLabel.Visible = 'on';
        %         end

    end

end

% linkaxes(ax(:));
% same_caxis(ax);
% title(t, currentDATA + "/" + parc_name{parc} + '/' + num2str(THR_LEVELS));
scfw(1300);
% print("v4_3a_" + currentDATA, '-dpng');

%% Fig S2b: Hub locations as a function of group reconstruction and density

verts = lh_inflated_verts;
faces = lh_faces;
currentDATA = ["ADJGroupWei", "ADJGroupDen", "ADJGroupConDenAltered", "ADJGroupDst"];
currentDATALabels = ["Weight", "CV", "Consistency", "Distance Dependence"];
THR_LEVELS = 1:11;
parc = 5;
rois = Scha_parcs.lh_scha200;
pipe = 7;

% ----

clear ax;
ax = gobjects(length(THR_LEVELS), length(currentDATA));



f = figure('Position', [500, 450, 160*length(currentDATA) , 1100], 'DefaultAxesFontSize', 30);
t1 = tiledlayout(1, length(currentDATA), "TileSpacing", 'none', 'Padding', 'none');



for jj = 1:length(currentDATA)

    DATA = eval(currentDATA(jj));


    %     currentSA = mean(PARC_SA{parc(jj)}, 1);


    t2 = tiledlayout(t1, length(THR_LEVELS), 1, 'TileSpacing', 'none', 'Padding', 'compact');
    t2.Layout.Tile = jj;

    for ii = 1:11


        current = DATA{ORDERED_INDS{parc}(pipe), THR_LEVELS(ii)};
        w = sum(current);
        [~, worder] = sort(w, 'ascend');

        ax(ii,jj) = nexttile(t2); axis(ax(ii,jj), 'off');

        t3 = tiledlayout(t2, 1, 2, 'TileSpacing', 'none', 'Padding', 'compact');
        t3.Layout.Tile = ii;


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

        %         if ii == 1
        %             ylabel(ax(ii,jj), parc_name{parc(jj)});
        %             ax(ii,jj).YLabel.Visible = 'on';
        %         end

    end

end

% linkaxes(ax(:));
% same_caxis(ax);
% title(t, currentDATA + "/" + parc_name{parc} + '/' + num2str(THR_LEVELS));
scfw(1200);
% print("v4_3a_" + currentDATA, '-dpng');