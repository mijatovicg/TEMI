function TE_rate = function_TE_rate(Cx, Jx, Cu, Ju, lambda_av_x, k_global)

% This function estimates the transfer entropy rate

% Input parameters:
% target history embeddings built at target events: Cx
% target history embeddings built at random events: Cu 
% joint history embeddings built at target events: Jx
% joint history embeddings built at random events: Ju
% average firing rate of a target spike train: lambda_av_x
% minimum number of nearest neighbors to consider it any search space: k_global

% Output parameter: 
% Transfer entropy rate: TE_rate

%% input parsing & validation
if ~exist('metric','var') % metric
    metric = 'maximum';
end

l_param = size(Cx, 2); % number of columns corresponds to length of embedding

%% application of k-nearest neighbor estimator
% Choose the maximum distance dd_max_C in lower dimension space and dd_max_J in higher dimension space

% sintax: [index, distance] = nn_search(pointset, atria, query_indices, k, exclude)
% sintax: [index, distance] = nn_search(pointset, atria, query_points, k)

atria_Cx = nn_prepare(Cx, metric); % internal search
[~, distances_Cxx] = nn_search(Cx, atria_Cx, (1:size(Cx,1))', k_global, 0); % sintax with query indices to exclude self-match since internal search is present!
dd_Cx = distances_Cxx(:, k_global); % % take distance to k_global_th neihgbor

atria_Cu = nn_prepare(Cu, metric); % external search
[~, distances_Cuu] = nn_search(Cu, atria_Cu, Cx, k_global); % sintax with query points
dd_Cu = distances_Cuu(:, k_global); % take distance to k_global_th neihgbor

%==========================================================================
dd_max_C = max(dd_Cx, dd_Cu); % maximum distances in lower dimension space
%==========================================================================

atria_Jx = nn_prepare(Jx, metric); % internal search
[~, distances_Jxx] = nn_search(Jx, atria_Jx, (1:size(Jx,1))', k_global, 0); 
dd_Jx = distances_Jxx(:, k_global);

atria_Ju = nn_prepare(Ju, metric); % external search
[~, distances_Juu] = nn_search(Ju, atria_Ju, Jx, k_global);
dd_Ju = distances_Juu(:, k_global); 

%==========================================================================
dd_max_J = max(dd_Jx, dd_Ju); % maximum distances in higher dimension space
%==========================================================================

%% range search, take different! number of neighbors for each space
% sintax: [count, neighbors] = range_search(pointset, atria, query_indices, r, exclude)
% sintax: [count, neighbors] = range_search(pointset, atria, query_points, r)

%--------------------------------------------------------------------------
% Cx 
[count_Cx, tmp_Cx] = range_search(Cx, atria_Cx, (1:size(Cx,1))', dd_max_C, 0); % internal search (Cx, Cx), sintax to exclude self-match!
%--------------------------------------------------------------------------
% Cu
if l_param == 1 % only if l_param == 1 add zeros column
    Cu_n = [Cu zeros(size(Cu,1),1)];
    Cx_n = [Cx zeros(size(Cx,1),1)];
    Cu = Cu_n;
    Cx = Cx_n;
end

[count_Cu, tmp_Cu] = range_search(Cu, atria_Cu, Cx, dd_max_C); % external search (Cx, Cu)
%--------------------------------------------------------------------------
% Jx 
[count_Jx, tmp_Jx] = range_search(Jx, atria_Jx, (1:size(Jx,1))', dd_max_J, 0); % sintax to exclude self-match!
%--------------------------------------------------------------------------
% Ju 
[count_Ju, tmp_Ju] = range_search(Ju, atria_Ju, Jx, dd_max_J);

%% with the different number of neighbors (count_Cx, count_Cu, count_Jx and count_Ju) search for the distances!
%--------------------------------------------------------------------------
distances_Cx = []; % internal search
for ss = 1 : size(Cx,1)
    [~, tmp] = nn_search(Cx, atria_Cx, ss, count_Cx(ss), 0);
    distances_Cx = [distances_Cx; tmp(count_Cx(ss))]; % take count_Cx(ss)_th distance!
end
%--------------------------------------------------------------------------
distances_Jx = []; % internal search
for ss = 1 : size(Jx,1)
    [~, tmp] = nn_search(Jx, atria_Jx, ss, count_Jx(ss), 0);
    distances_Jx = [distances_Jx; tmp(count_Jx(ss))]; % take count_Jx(ss)_th distance!
end
%--------------------------------------------------------------------------
distances_Cu = []; % external search
for ss = 1 : size(Cx,1)
    [~, tmp] = nn_search(Cu, atria_Cu, Cx(ss, :), count_Cu(ss));
    distances_Cu = [distances_Cu; tmp(count_Cu(ss))]; % take count_Cy(ss)_th distance!
end
%--------------------------------------------------------------------------
distances_Ju = []; % external search
for ss = 1 : size(Jx,1)
    [~, tmp] = nn_search(Ju, atria_Ju, Jx(ss, :), count_Ju(ss));
    distances_Ju = [distances_Ju; tmp(count_Ju(ss))]; % take count_Ju(ss)_th distance!
end

%% take twice distance
distances_Cx = 2*distances_Cx;
distances_Cu = 2*distances_Cu;
distances_Jx = 2*distances_Jx;
distances_Ju = 2*distances_Ju;

%% transfer entropy rate estimation
TE = mean(l_param*(log(distances_Cx) - log(distances_Cu)) - 2*l_param*(log(distances_Jx) - log(distances_Ju)) - psi(count_Cx) + psi(count_Cu) + psi(count_Jx) - psi(count_Ju));
TE_rate = TE * lambda_av_x;

