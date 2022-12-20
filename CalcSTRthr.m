function STRDATA = CalcSTRthr(FILEDIR,save_output_location,IGNORE_SUBCORTEX)

if nargin < 3
    IGNORE_SUBCORTEX = 0;
end
% This function calculates the correlation between edge weight and some
% measure of motion across participants.
%
% INPUTS:
% FILEDIR = location of where the files to load in are located. Assumes
% they have the format "Pipeline_$.mat' where $ is the pipeline number
%
% save_output_location = set the location of where output is saved
% (defaults to current directory)
%
% IGNORE_SUBCORTEX = set to 1 to remove the subcortex from strength calculations
%
% OUTPUTS:
%
% STRDATA = a structure of all data related to node degree and strength.
% See below for explanation of all fields in this structure. Each contains
% a cell/vector of some property for each pipeline.


% Apply a basic consistency threshold to all connectoms before any other
% threshold is applied. This threshold will remove any edges which have a
% consistency < edgeconthr. Probably want to leave it at 0
edgeconthr = 0;

if nargin < 1
    save_output_location = [];
end

% Get information on the COMBINATIONS
load('COMBINATIONS_MATRIX.mat')

% Get the mean Euclidean diststances between ROIs across participants
load('SApipes_meandist.mat','MeanDists')

NPipes = size(COMBINATIONS,1);

% Thresholds to check when calculating strength
threshs = [.05 .075 .1 .125 .15 .175 .2 .225 .25 .275 .3];

threshs_con = [1 .95 .90 .85 .80 .75 .7 .65 .6 .55 .5];

threshs_bins = [10 20 30 40 50 60 70 80 90 100 110];

Nthr = length(threshs);

% A cell matrix where each cell is the node strength for a pipeline (rows)
% under a given edge based consistency threshold (columns). This keeps
% edges which have a minimum consistency (present across subjects)
STRDATA.STRcon = cell(NPipes,Nthr);

% A cell matrix where each cell is the node degree for a pipeline (rows)
% under a given edge based consistency threshold (columns). This keeps
% edges which have a minimum consistency (present across subjects)
STRDATA.DEGcon = cell(NPipes,Nthr);

% A cell matrix where each cell is the node strength for a pipeline (rows)
% under a given edge based variability threshold (columns). Keeps a
% proportion of all edges with the smallest coefficient of variation
STRDATA.STRvar = cell(NPipes,Nthr);

% A cell matrix where each cell is the node degree for a pipeline (rows)
% under a given edge based variability threshold (columns). Keeps a
% proportion of all edges with the smallest coefficient of variation
STRDATA.DEGvar = cell(NPipes,Nthr);

% A cell matrix where each cell is the node strength for a pipeline (rows)
% under a given edge based variability threshold (columns). Keeps the edges
% with the smallest coefficient of variation up to a given density
STRDATA.STRden = cell(NPipes,Nthr);

% A cell matrix where each cell is the node degree for a pipeline (rows)
% under a given edge based variability threshold (columns). Keeps the edges
% with the smallest coefficient of variation up to a given density
STRDATA.DEGden = cell(NPipes,Nthr);

% A cell matrix where each cell is the node strength for a pipeline (rows)
% under a given edge based variability threshold (columns). Keeps the edges
% with the largest weight up to a given density
STRDATA.STRwei = cell(NPipes,Nthr);

% A cell matrix where each cell is the node degree for a pipeline (rows)
% under a given edge based variability threshold (columns). Keeps the edges
% with the largest weight up to a given density
STRDATA.DEGwei = cell(NPipes,Nthr);

% A cell matrix where each cell is the node strength for a pipeline (rows)
% under a given edge based variability threshold (columns). Applies a
% threshold across different edge bin lengths
STRDATA.STRdst = cell(NPipes,Nthr);

% A cell matrix where each cell is the node degree for a pipeline (rows)
% under a given edge based variability threshold (columns). Applies a
% threshold across different edge bin lengths
STRDATA.DEGdst = cell(NPipes,Nthr);

% A cell matrix where each cell is the node strength for a pipeline (rows)
% under a given edge based variability threshold (columns)
STRDATA.STRconden = cell(NPipes,Nthr);

% A cell matrix where each cell is the node degree for a pipeline (rows)
% under a given edge based variability threshold (columns)
STRDATA.DEGconden = cell(NPipes,Nthr);

STRDATA.threshs = threshs;

STRDATA.threshs_con = threshs_con;

STRDATA.threshs_bins = threshs_bins;

for ITER = 1:size(COMBINATIONS,1)

    %% Extract data
tic
clear adjs lens con

load([FILEDIR,'Pipeline_',num2str(ITER),'.mat'],'adjs','lens')

% We need to know the number of regions in each parcellation. However
% sometimes we want to ignore the subcortex. So the following loop
% determines which parcellation is being used (based on the COMBINATIONS
% matrix) and indicates which regions to use (ROIS)

PARC_IND = COMBINATIONS(ITER,7);

if IGNORE_SUBCORTEX
    if PARC_IND == 2
        ROIS = [1:100 111:210];
    elseif PARC_IND == 1
        ROIS = [1:34 42:75];
    elseif PARC_IND == 3
        ROIS = [1:180 191:370];
    elseif PARC_IND == 4
        ROIS = [1:250 261:510];
    elseif PARC_IND == 5
        ROIS = [1:100 111:210];
    elseif PARC_IND == 6
        ROIS = [1:250 261:510];
    end
    % Get the connectivity matrix for the ROIs
    for i = 1:length(adjs)
        adjs{i} = adjs{i}(ROIS,ROIS);
        lens{i} = lens{i}(ROIS,ROIS);
    end

else
    ROIS = 1:length(adjs{1});
end

% Calculate the consistency of each edge across subjects (proporartion of
% times that edge exists)
[~,~,~,con] = connectomeGroupThreshold(adjs,0);

Nsubs = length(adjs);
NNodes = length(con);

% Make a group binary mask
groupMask = con;
groupMask(groupMask<edgeconthr) = 0;
groupMask(groupMask~=0) = 1;

% Get the upper triangle
groupMaskTriu = triu(groupMask,1);

% Get an index for each unique edge
[indEdge] = find(groupMaskTriu(:)>0);

NEdges = length(indEdge);

% Prepare matrices for a subject by edge matrix
subjEdges = zeros(Nsubs, NEdges);
subjEdgesLength = zeros(Nsubs, NEdges);

% Ws is a 3D matrix for the subjects connectomes
Ws = zeros(NNodes,NNodes,Nsubs);
% den is subject densities
den = zeros(Nsubs,1);
% total_STR is the sum of all edge weights
total_STR = zeros(Nsubs,1);
% means and sds for each subjects edge weights
EdgeWeight_mean = zeros(Nsubs,1);
EdgeWeight_sd = zeros(Nsubs,1);

% loop over subjects
for subj=1:Nsubs

    % Get a vector of the subjects edge weights and lengths
    adjVect = adjs{subj}(:);
    lengthVect = lens{subj}(:);

    % Pull out the weights for all the edges identified by indEdge
    edge_weights = adjVect(indEdge);

    % Insert edges into edge matrix
    subjEdges(subj,:) = edge_weights;

    % Insert lengths into edge length matrix
    subjEdgesLength(subj,:) = lengthVect(indEdge);

    % Apply the binary mask
    W = adjs{subj}.*groupMask;

    % Insert into matrix
    Ws(:,:,subj) = W;

    % Calculate density
    Nnodes = size(W,1);
    Nedges = nnz(W(~isnan(W)))/2;
    den(subj,1) = Nedges/((Nnodes^2-Nnodes)/2);

    % Calculate total strength, mean and sd of edge weights
    total_STR(subj,1) = nansum(edge_weights);
    EdgeWeight_mean(subj,1) = nanmean(edge_weights(edge_weights>0));
    EdgeWeight_sd(subj,1) = nanstd(edge_weights(edge_weights>0));

end

% Calculate the coefficient of variation (CV) for each edge
[~,Wcv] = threshold_consistency(Ws, 1);
% Get the upper triangle and index of each unique edge
[Wcv_vec,Wcv_ind] = triu2vec(Wcv,1);
% Sort the edges from smallest to largest CV
[~,Wcv_ind_order] = sort(Wcv_vec,'ascend');

% Apply different thresholds at different levels
for thr = 1:length(threshs)
    % Apply a consistency threshold (keep edges which are present in at
    % least threshs(thr) subjects
    [~, ~, AdjMaskCon] = connectomeGroupThreshold(adjs,threshs_con(thr));

    % Apply a threshold based on CV (keep the edges with the smallest CV),
    % will keep a proportion of all existing connections
    AdjMaskVar = threshold_consistency(Ws, threshs(thr));

    % Apply a threshold based on CV, will threshold to a certain density
    % (by only keeping those with the smallest CV)
    Den_Wcv_inds = Wcv_ind(Wcv_ind_order(1:floor(length(Wcv_ind)*threshs(thr))));
    AdjMaskDen = zeros(NNodes);
    AdjMaskDen(Den_Wcv_inds) = 1;
    AdjMaskDen = AdjMaskDen + AdjMaskDen';

    % Do length based thresholding. See fcn_group_bins for more
    % information.
    DISTS = MeanDists{COMBINATIONS(ITER,7)};
    DISTS = DISTS(ROIS,ROIS);
    hemiid = [ones(length(DISTS)/2,1); ones(length(DISTS)/2,1)*2];
    ADJGroupDst = double(fcn_group_bins(Ws,DISTS,hemiid,threshs_bins(thr))>0);

    % Apply a weight based threshold by keeping the strongest edges up to a
    % certain density for each subject
    WsThr = zeros(NNodes,NNodes,Nsubs);
    for i = 1:Nsubs
        w = adjs{subj};
        [w_vec,w_ind] = triu2vec(w,1);
        [~,w_ind_order] = sort(w_vec,'descend');
        w_thr_ind = w_ind(w_ind_order(1:floor(length(w_ind)*threshs(thr))));
        w_thr = zeros(NNodes);
        w_thr(w_thr_ind) = 1;
        WsThr(:,:,i) = (w_thr + w_thr').*w;
    end

    % Get the average edge weight across subjects (unthresholded)
    Ws_mean = nanmean(Ws,3);

    % Keeps the most consistent edges to reach a given threshold
    % Sorts by strength if consistency is tied
    [w_vec,w_ind] = triu2vec(con,1);
    s_vec = triu2vec(Ws_mean,1);
    [~,w_ind_order] = sortrows([w_vec, s_vec],'descend');
    w_thr_ind = w_ind(w_ind_order(1:floor(length(w_ind)*threshs(thr))));
    w_thr = zeros(NNodes);
    w_thr(w_thr_ind) = 1;
    AdjMaskConDen = (w_thr + w_thr');
    STRDATA.ADJGroupConDen{ITER,thr} = AdjMaskConDen.*Ws_mean;
    ConDenAdjs = zeros(NNodes,NNodes,Nsubs);

    for i = 1:Nsubs
        ConDenAdjs(:,:,i) = AdjMaskConDen.*adjs{i};
    end

    % Applies the binary thresholded mask to the mean connectome
STRDATA.ADJGroupCon{ITER,thr} = Ws_mean.*AdjMaskCon;
STRDATA.ADJGroupVar{ITER,thr} = Ws_mean.*AdjMaskVar;
STRDATA.ADJGroupDen{ITER,thr} = Ws_mean.*AdjMaskDen;
STRDATA.ADJGroupDst{ITER,thr} = Ws_mean.*ADJGroupDst;

% Applies a density threshold to the mean connectome to keep the strongest
% edges up to a given density.

% For the weighted thresholding, the group connectome isn't created by
% averaging the subject weighted thresholded matrices, but by appling the
% weight threshold to the mean weighted matrix. This is different to the
% other thresholding options as they defining a binary matrix which is
% applied to the mean matrix. This might be wrong then...
[w_vec,w_ind] = triu2vec(nanmean(Ws,3),1);
[~,w_ind_order] = sort(w_vec,'descend');
w_thr_ind = w_ind(w_ind_order(1:floor(length(w_ind)*threshs(thr))));
w_thr = zeros(NNodes);
w_thr(w_thr_ind) = 1;
STRDATA.ADJGroupWei{ITER,thr} = (w_thr + w_thr').*Ws_mean;

% Apply the group mask to each individuals connectome
ConAdjs = Ws.*repmat(AdjMaskCon,1,1,Nsubs);
VarAdjs = Ws.*repmat(AdjMaskVar,1,1,Nsubs);
DenAdjs = Ws.*repmat(AdjMaskDen,1,1,Nsubs);
DstAdjs = Ws.*repmat(ADJGroupDst,1,1,Nsubs);

WeiAdjs = WsThr;

% The variables are doubled up instead of directly being assigned to a cell
% because it made troubleshooting easier

Str = squeeze(nansum(ConAdjs,2))';
Deg = squeeze(nansum(ConAdjs>0,2))';

Str2 = squeeze(nansum(VarAdjs,2))';
Deg2 = squeeze(nansum(VarAdjs>0,2))';

Str3 = squeeze(nansum(DenAdjs,2))';
Deg3 = squeeze(nansum(DenAdjs>0,2))';

Str4 = squeeze(nansum(WeiAdjs,2))';
Deg4 = squeeze(nansum(WeiAdjs>0,2))';

Str5 = squeeze(nansum(DstAdjs,2))';
Deg5 = squeeze(nansum(DstAdjs>0,2))';

Str6 = squeeze(nansum(ConDenAdjs,2))';
Deg6 = squeeze(nansum(ConDenAdjs>0,2))';

STRDATA.STRcon{ITER,thr} = Str;
STRDATA.DEGcon{ITER,thr} = Deg;
STRDATA.STRvar{ITER,thr} = Str2;
STRDATA.DEGvar{ITER,thr} = Deg2;

STRDATA.STRden{ITER,thr} = Str3;
STRDATA.DEGden{ITER,thr} = Deg3;

STRDATA.STRwei{ITER,thr} = Str4;
STRDATA.DEGwei{ITER,thr} = Deg4;

STRDATA.STRdst{ITER,thr} = Str5;
STRDATA.DEGdst{ITER,thr} = Deg5;

STRDATA.STRconden{ITER,thr} = Str6;
STRDATA.DEGconden{ITER,thr} = Deg6;

end

timetaken = toc;
fprintf('Completed %d/60 in %.4f seconds\n',ITER,timetaken)

end

STRDATA.COMBINATIONS = COMBINATIONS;

if ~isempty(save_output_location)

   save([save_output_location],'-struct','STRDATA','-v7.3')

end


function [merged_struct] = MergeStructs(struct_a,struct_b)
%%if one of the structres is empty do not merge
if isempty(struct_a)
    merged_struct=struct_b;
    return
end
if isempty(struct_b)
    merged_struct=struct_a;
    return
end
%%insert struct a
merged_struct=struct_a;
%%insert struct b
size_a=length(merged_struct);
for j=1:length(struct_b)
    f = fieldnames(struct_b);
    for i = 1:length(f)
        merged_struct(size_a+j).(f{i}) = struct_b(j).(f{i});
    end
end
