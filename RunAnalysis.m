% This is what is used to make suammry figures for each pipeline.
% Heatmaps, strength/degree vs SA, etc.
% I'll do my best to explain each section. WARNING: long

addpath(genpath('./'))

% Load in surface area for each region in each parcellation (no
% subcortical)
load('Parc_SA.mat')

% Load in Euclidean distances between ROIs for each parcellation
load('SApipes_meandist.mat')

% Load in the COMBINATIONS matrix. This indicates for each pipeline (row)
% what processing option was specified.
load('COMBINATIONS_MATRIX.mat')

% This will make the indexing information you need to pull out relevant
% data for each pipeline

for parc = 1:6
    [ORDERED_INDS{parc},ORDERED_MATRIX{parc},LABELS{parc}] = FindPipelineCombinations([1 0 0 0 1 1 parc],[7 1 6 2 4 3 5],1);
end

% More explicitly define the names of the processing steps
pipeline_titles = {'ACT/dynamic/FACT','GWM/dynamic/FACT','ACT/WM/FACT','GWM/WM/FACT','ACT/GMWMI/FACT',...
    'ACT/dynamic/iFOD2','GWM/dynamic/iFOD2','ACT/WM/iFOD2','GWM/WM/iFOD2','ACT/GMWMI/iFOD2'};



%% --------------------------------------------------------------------- %%

% Calculate the strength/degree across all processing pipelines. These will
% take a long time to run FYI. If they are already generated, comment this part out!
% CalcSTRthr('./DATA','./STRDATA_SUBCORT.mat',0);
% CalcSTRthr('./DATA','./STRDATA_NOSUBCORT.mat',1);



%% --------------------------------------------------------------------- %%


% Can do both with and without cortex by setting DATATYPE to 1 (with subcortex)
% or 2 (without subcortex)

% Select which threshold type and whether to use degree or strengh to use (1:12)
% Running over all of them may take a while depending on how many threshold levels you have selected
% 1 = DEGden
% 2 = STRden
% 3 = DEGcon
% 4 = STRcon
% 5 = DEGvar
% 6 = STRvar
% 7 = DEGwei
% 8 = STRwei
% 9 = DEGdst
% 10 = STRdst
% 11 = DEGconden
% 12 = STRconden
THR_TYPES_2_USE = [1 2];

% Select the threshold levels to use (1:11). The con threshold type uses thresholds corresponding to that in
% thr_strings_consistency, all others apart from dst use those in thr_strings_density, and dst uses those in
% thr_strings_dst

THR_LEVELS = [1 5 7 11];

thr_strings_consistency = {'100','95','90','85','80','75','70','65','60','55','50'};
thr_strings_density = {'5','7.5','10','12.5','15','17.5','20','22.5','25','27.5','30'};
thr_strings_dst = {'10','20','30','40','50','60','70','80','90','100','110'};


%% --------------------------------------------------------------------- %%

for DATATYPE = 1:2

    if DATATYPE == 1
        SAVELOC = './STRDATA_SUBCORT';
        % Probably should load these in as structures. Neater that way
        load('STRDATA_SUBCORT.mat')
        nosubcortex = 0;

    elseif DATATYPE == 2
        SAVELOC = './STRDATA_NOSUBCORT';
        load('STRDATA_NOSUBCORT.mat')
        nosubcortex = 1;
    end

    % Make the directory with which to save the many figures too
    mkdir(SAVELOC)

    % Get the name of the parcellation
    parc_name = Parcellation;

    % Loop over parcellations
    for parc = 1:6
        % Depending on if the subcortex is being included or not, get the
        % number of nodes, and which regions are cortical
        if parc == 2
            n = 220;
            CORT = [1:100 111:210];
            RANGE = 1:220;
        elseif parc == 1
            n = 82;
            CORT = [1:34 42:75];
            RANGE = 1:82;
        elseif parc == 3
            n = 380;
            CORT = [1:180 191:370];
            RANGE = 1:380;
        elseif parc == 4
            n = 520;
            CORT = [1:250 261:510];
            RANGE = 1:520;
        elseif parc == 5
            n = 220;
            CORT = [1:100 111:210];
            RANGE = [1:100 201:310 411:420];
        elseif parc == 6
            n = 520;
            CORT = [1:250 261:510];
            RANGE = [1:250 351:610 711:720];
        end

        if nosubcortex == 1
            if parc == 1
                CORT = 1:(n-14);
            else
                CORT = 1:(n-20);
            end
        end

        % Loop over each pipeline
        for i = 1:length(ORDERED_INDS{parc})
            % Loop over each trehsold level
            for j = 1:length(threshs)
                % This goes through each pipeline, for each threshold level and
                % calculates the skewness of the degree/strength distribution
                STR = STRcon{ORDERED_INDS{parc}(i),j};
                DEG = DEGcon{ORDERED_INDS{parc}(i),j};
                skew = skewness(STR');
                skew_deg = skewness(DEG');

                indiv_skew_str_con(j,i) = mean(skew);
                indiv_skew_deg_con(j,i) = mean(skew_deg);

                STR = STRvar{ORDERED_INDS{parc}(i),j};
                DEG = DEGvar{ORDERED_INDS{parc}(i),j};
                skew = skewness(STR');
                skew_deg = skewness(DEG');

                indiv_skew_str_var(j,i) = mean(skew);
                indiv_skew_deg_var(j,i) = mean(skew_deg);


                STR = STRden{ORDERED_INDS{parc}(i),j};
                DEG = DEGden{ORDERED_INDS{parc}(i),j};
                skew = skewness(STR');
                skew_deg = skewness(DEG');

                indiv_skew_str_den(j,i) = mean(skew);
                indiv_skew_deg_den(j,i) = mean(skew_deg);

                STR = STRwei{ORDERED_INDS{parc}(i),j};
                DEG = DEGwei{ORDERED_INDS{parc}(i),j};
                skew = skewness(STR');
                skew_deg = skewness(DEG');

                indiv_skew_str_wei(j,i) = mean(skew);
                indiv_skew_deg_wei(j,i) = mean(skew_deg);

                STR = STRdst{ORDERED_INDS{parc}(i),j};
                DEG = DEGdst{ORDERED_INDS{parc}(i),j};
                skew = skewness(STR');
                skew_deg = skewness(DEG');

                indiv_skew_str_dst(j,i) = mean(skew);
                indiv_skew_deg_dst(j,i) = mean(skew_deg);

                STR = STRconden{ORDERED_INDS{parc}(i),j};
                DEG = DEGconden{ORDERED_INDS{parc}(i),j};
                skew = skewness(STR');
                skew_deg = skewness(DEG');

                indiv_skew_str_conden(j,i) = mean(skew);
                indiv_skew_deg_conden(j,i) = mean(skew_deg);

            end
        end

        % Once we have calculated all the skewness values, print a figure for
        % them!

        figure
        PlotMatrixWithLegend(indiv_skew_str_con,ORDERED_MATRIX{parc},LABELS{parc},'Skewness',threshs_con)
        print([SAVELOC,'/Skewness_',num2str(parc_name{parc}),'_str_consistency.png'],'-dpng')

        figure
        PlotMatrixWithLegend(indiv_skew_deg_con,ORDERED_MATRIX{parc},LABELS{parc},'Skewness',threshs_con)
        print([SAVELOC,'/Skewness_',num2str(parc_name{parc}),'_deg_consistency.png'],'-dpng')

        figure
        PlotMatrixWithLegend(indiv_skew_str_var,ORDERED_MATRIX{parc},LABELS{parc},'Skewness',threshs)
        print([SAVELOC,'/Skewness_',num2str(parc_name{parc}),'_str_variance.png'],'-dpng')

        figure
        PlotMatrixWithLegend(indiv_skew_deg_var,ORDERED_MATRIX{parc},LABELS{parc},'Skewness',threshs)
        print([SAVELOC,'/Skewness_',num2str(parc_name{parc}),'_deg_variance.png'],'-dpng')

        figure
        PlotMatrixWithLegend(indiv_skew_str_den,ORDERED_MATRIX{parc},LABELS{parc},'Skewness',threshs)
        print([SAVELOC,'/Skewness_',num2str(parc_name{parc}),'_str_density.png'],'-dpng')

        figure
        PlotMatrixWithLegend(indiv_skew_deg_den,ORDERED_MATRIX{parc},LABELS{parc},'Skewness',threshs)
        print([SAVELOC,'/Skewness_',num2str(parc_name{parc}),'_deg_density.png'],'-dpng')

        figure
        PlotMatrixWithLegend(indiv_skew_str_wei,ORDERED_MATRIX{parc},LABELS{parc},'Skewness',threshs)
        print([SAVELOC,'/Skewness_',num2str(parc_name{parc}),'_str_weight.png'],'-dpng')

        figure
        PlotMatrixWithLegend(indiv_skew_deg_wei,ORDERED_MATRIX{parc},LABELS{parc},'Skewness',threshs)
        print([SAVELOC,'/Skewness_',num2str(parc_name{parc}),'_deg_weight.png'],'-dpng')

        figure
        PlotMatrixWithLegend(indiv_skew_str_dst,ORDERED_MATRIX{parc},LABELS{parc},'Skewness',threshs_bins)
        print([SAVELOC,'/Skewness_',num2str(parc_name{parc}),'_str_dst.png'],'-dpng')

        figure
        PlotMatrixWithLegend(indiv_skew_deg_dst,ORDERED_MATRIX{parc},LABELS{parc},'Skewness',threshs_bins)
        print([SAVELOC,'/Skewness_',num2str(parc_name{parc}),'_deg_dst.png'],'-dpng')

        figure
        PlotMatrixWithLegend(indiv_skew_str_conden,ORDERED_MATRIX{parc},LABELS{parc},'Skewness',threshs_bins)
        print([SAVELOC,'/Skewness_',num2str(parc_name{parc}),'_str_conden.png'],'-dpng')

        figure
        PlotMatrixWithLegend(indiv_skew_deg_conden,ORDERED_MATRIX{parc},LABELS{parc},'Skewness',threshs_bins)
        print([SAVELOC,'/Skewness_',num2str(parc_name{parc}),'_deg_conden.png'],'-dpng')
    end

    % This loops over different thresholding options and makes a number of
    % matrices. Does degree and strength distributions seperately
    for T = THR_TYPES_2_USE
        if T == 1
            DATA = DEGden;
            DATAname = 'Degree';
            savename = 'DEGDEN';
        elseif T == 2
            DATA = STRden;
            DATAname = 'Strength';
            savename = 'STRDEN';
        elseif T == 3
            DATA = DEGcon;
            DATAname = 'Degree';
            savename = 'DEGCON';
        elseif T == 4
            DATA = STRcon;
            DATAname = 'Strength';
            savename = 'STRCON';
        elseif T == 5
            DATA = DEGvar;
            DATAname = 'Degree';
            savename = 'DEGVAR';
        elseif T == 6
            DATA = STRvar;
            DATAname = 'Strength';
            savename = 'STRVAR';
        elseif T == 7
            DATA = DEGwei;
            DATAname = 'Degree';
            savename = 'DEGWEI';
        elseif T == 8
            DATA = STRwei;
            DATAname = 'Strength';
            savename = 'STRWEI';
        elseif T == 9
            DATA = DEGdst;
            DATAname = 'Degree';
            savename = 'DEGDST';
        elseif T == 10
            DATA = STRdst;
            DATAname = 'Strength';
            savename = 'STRDST';
        elseif T == 11
            DATA = DEGconden;
            DATAname = 'Degree';
            savename = 'DEGCONDEN';
        elseif T == 12
            DATA = STRconden;
            DATAname = 'Strength';
            savename = 'STRCONDEN';
        end

        for thr = THR_LEVELS
            if ismember(T,[3 4])

                thr_string = thr_strings_consistency{thr};

            elseif ismember(T,[9 10])
                thr_string = thr_strings_dst{thr};
            else
                thr_string = thr_strings_density{thr};
            end

            for parc = 1:6
                figure('Position',[933 885 1761 603])
                if parc == 2
                    n = 220;
                    CORT = [1:100 111:210];
                elseif parc == 1
                    n = 82;
                    CORT = [1:34 42:75];
                elseif parc == 3
                    n = 380;
                    CORT = [1:180 191:370];
                elseif parc == 4
                    n = 520;
                    CORT = [1:250 261:510];
                elseif parc == 5
                    n = 220;
                    CORT = [1:100 111:210];
                elseif parc == 6
                    n = 520;
                    CORT = [1:250 261:510];
                end

                % Indicate the total number of regions if the subcortex is included or not.
                % This overwrites the CORT variable defined above if required
                if nosubcortex == 1
                    if parc == 1
                        % In this parcellation there are 14 subcortical regions
                        n = n-14;

                    else
                        % In these parcellation there are 20 subcortical regions
                        n = n-20;
                    end
                    CORT = 1:n;
                end

                % Across individuals calculate the mean strength of each region and
                % plot the histogram of these mean strengths
                mean_strength_rank = zeros(n,10);
                var_strength = zeros(n,10);
                % Loop over pipelines
                for i = 1:length(ORDERED_INDS{parc})

                    DEG = DATA{ORDERED_INDS{parc}(i),thr};
                    mean_strength_rank(:,i) = mean(DEG,1);
                    var_strength(:,i) = var(DEG,1);

                    subplot(2,5,i)
                    histogram(mean_strength_rank(:,i))
                    title(pipeline_titles{i})
                    xlabel(['Mean ',DATAname])
                end
                print([SAVELOC,'/Histogram_',num2str(parc_name{parc}),'_',savename,'_',thr_string,'.png'],'-dpng')
                meanSA = mean(PARC_SA{parc});
                varSA = var(PARC_SA{parc});

                meanDist = MeanDists{parc}(CORT,CORT);

                % This calculates the correlation between mean strength and mean surface
                % area for each region. Also plots these against each other, colouring each
                % point by the mean distance of that region to others
                figure('Position',[933 885 1761 603])
                for i = 1:length(ORDERED_INDS{parc})
                    subplot(2,5,i)
                    scatter(mean_strength_rank(CORT,i),meanSA,20,mean(meanDist),'filled')
                    cp = corr(mean_strength_rank(CORT,i),meanSA');
                    cs = corr(mean_strength_rank(CORT,i),meanSA','Type','Spearman');
                    title({pipeline_titles{i},[' R = ',num2str(round(cp,4)),' Rho = ',num2str(round(cs,4))]})
                    hold on
                    % Fit a linear trend
                    p = polyfit(mean_strength_rank(CORT,i),meanSA,1);
                    f = polyval(p,mean_strength_rank(CORT,i));
                    plot(mean_strength_rank(CORT,i),f,'r')
                    xlabel(['Mean ',DATAname])
                    ylabel('Mean Surface Area')
                end

                c = colorbar('Position',[0.921307871259802,0.15257048092869,0.016227619938381,0.70978441127695]);
                c.FontSize = 16;
                c.Label.String = 'Mean distance';
                print([SAVELOC,'/MeanSA_',num2str(parc_name{parc}),'_',savename,'_',thr_string,'.png'],'-dpng')


                % Same as above but for ROI strength and surface area variability
                figure('Position',[933 885 1761 603])

                for i = 1:length(ORDERED_INDS{parc})
                    subplot(2,5,i)
                    scatter(var_strength(CORT,i),meanSA,20,mean(meanDist),'filled')
                    cp = corr(var_strength(CORT,i),meanSA');
                    cs = corr(var_strength(CORT,i),meanSA','Type','Spearman');
                    title({pipeline_titles{i},[' R = ',num2str(round(cp,4)),' Rho = ',num2str(round(cs,4))]})
                    hold on
                    p = polyfit(var_strength(CORT,i),meanSA,1);
                    f = polyval(p,var_strength(CORT,i));
                    plot(var_strength(CORT,i),f,'r')
                    xlabel(['Variance ',DATAname])
                    ylabel('Variance Surface Area')
                end

                c = colorbar('Position',[0.921307871259802,0.15257048092869,0.016227619938381,0.70978441127695]);
                c.FontSize = 16;
                c.Label.String = 'Mean distance';
                print([SAVELOC,'/VarSA_',num2str(parc_name{parc}),'_',savename,'_',thr_string,'.png'],'-dpng')

                % Makes a correlation matrix of the similarity in degree/strength distribution
                % between pipelines, clusters similar pipelines together

                CSTR = corr(mean_strength_rank,'Type','Spearman');
                RunClusterPipelineProp(CSTR,2,1,ORDERED_MATRIX{parc},LABELS{parc},[DATAname,' correlation']);
                print([SAVELOC,'/CorrMat_',num2str(parc_name{parc}),'_',savename,'_',thr_string,'.png'],'-dpng')

                % Does as above but for a partial correlation controlling for mean surface
                % area

                CSTR = partialcorr(mean_strength_rank(CORT,:),meanSA','Type','Spearman');
                RunClusterPipelineProp(CSTR,2,1,ORDERED_MATRIX{parc},LABELS{parc},[DATAname,' partial correlation']);
                print([SAVELOC,'/PartialCorrMat_',num2str(parc_name{parc}),'_',savename,'_',thr_string,'.png'],'-dpng')

                close all

            end

        end

    end

end
