function MI_rate = function_MI_rate(Cu_x, Ju_x, Cu_y, lambda_av_u, k_global)

% This function estimates the dynamic mutual information rate

% Input parameters:
% history embeddings built at random events of one spike train: Cu_x 
% history embeddings built at random events of another spike train: Cu_y
% joint history embeddings built at random events: Ju_x
% average firing rate of a stream of randomly generated events: lambda_av_u
% minimum number of nearest neighbors to consider it any search space: k_global

% Output parameter: 
% Dynamic mutual information rate:MI_rate


%% input parsing & validation
if ~exist('metric','var') % metric
    metric = 'maximum';
end

Nu = size(Cu_x, 1); % number of rows corresponds to the total number of random events

%% distances in higher dimensional space
Ju = Ju_x; % (same result if Ju_y is used)
atria_Ju = nn_prepare(Ju, metric); % internal search
[~, distances_Juu] = nn_search(Ju, atria_Ju, (1:size(Ju,1))', k_global, 0); % sintax with query indices to exclude self-match since internal search is present!
dd_Ju = distances_Juu(:, k_global); % % take distance to k_global_th neihgbor

%% number of neighbors in lower dimensional spaces
atria_Cu_x = nn_prepare(Cu_x, metric);
[count_Cu_x, tmp_Cu_x] = range_search(Cu_x, atria_Cu_x, (1:size(Cu_x,1))', dd_Ju, 0); % internal search

atria_Cu_y = nn_prepare(Cu_y, metric);
[count_Cu_y, tmp_Cu_y] = range_search(Cu_y, atria_Cu_y, (1:size(Cu_y,1))', dd_Ju, 0); % internal search

%% MI rate estimation
MI_rate = lambda_av_u*(psi(k_global) + log(Nu-1) - mean(psi(count_Cu_x) + psi(count_Cu_y)));

