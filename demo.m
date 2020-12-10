close all; clc;

% Here is demonstrated the framework for the transfer entropy (TE) rate and  
% dynamical mutual information (dMI) estimation, proposed in the paper 
% "An Information?Theoretic Framework to Measure the Dynamic Interaction between Neural Spike Trains" 
% by G. Mijatovic, Y. Antonacci, T. Loncar Turukalo, L. Minati and L. Faes,2020. 

% The framework presents an information-theoretic approach for the model-free, continuous-time estimation 
% of both undirected, symmetric (dMI) and directed, causal (TE) interactions between spike trains.

% Script demo.m simulates two coupled spike trains denoted as spike train 1
% (process X) and 2 (process Y) where unidirectional interaction X--?Y is 
% simulated setting tau=delta=0.1 s(see the simulation protocol in Sect. III in the main article).

% This script uses next functions: 
% function_embedding_vectors.m; 
% function_TE_rate; 
% function_MI_rate.m; 
% spiSeMe_surrogate_jodi.m 
% as a part of package 'TEMI'


%% input parameters
T = 300; % duration of spike trains
rate = 1; % average firing rate of both trains

l_param = 1; % length of embedding
Nx_Nu_ratio = 1; % ration between target and randomly generated events
k_global = 5; % minimum number of the nearest neighbors to consider it any set.

delta = 0.1; % paramater of simulation of coupled trains
Tau = delta; % paramater of simulation of coupled trains

n_surr = 100; % number of surrogates

%% generate coupled Poisson spike trains
spike_train1 = generatePoisson(rate, T);
U = 2*delta.*rand(numel(spike_train1),1) - delta; % jitter uniformly distribed in the range [-delta, delta]
spike_train2 = spike_train1 + Tau + U;

%% generate random events
N_spikes = max(numel(spike_train1), numel(spike_train2));
if numel(spike_train1) ==  N_spikes
    time_limit = max(spike_train1);
elseif numel(spike_train2) ==  N_spikes
    time_limit = max(spike_train2);
else
end
Nu = round(N_spikes*Nx_Nu_ratio); % ratio between number of spikes
random_events = sort(round(time_limit).*rand(Nu, 1)); % gererate radnom time axis by uniform distribution

%% TE rate 2 --> 1 (spike_train1: target, spike_train2: driver)
spike_train_target = spike_train1;
spike_train_driver = spike_train2;

tau  = spike_train_target(end, :) - spike_train_target(1, :); % total recording duration!
lambda_av_2_1 = numel(spike_train_target)/tau; % average firing rate
[Cx, Jx, Cu_x, Ju_x] = function_embedding_vectors(spike_train_target, spike_train_driver, l_param, random_events); % create embeddings
TE_rate_2_1 = function_TE_rate(Cx, Jx, Cu_x, Ju_x, lambda_av_2_1, k_global);

%% TE rate 1 --> 2 (spike_train2: target, spike_train1: driver)
spike_train_target = spike_train2;
spike_train_driver = spike_train1;

tau  = spike_train_target(end, :) - spike_train_target(1, :); % last spike, total recording duration!
lambda_av_1_2 = numel(spike_train_target)/tau; % average firing rate
[Cy, Jy, Cu_y, Ju_y] = function_embedding_vectors(spike_train_target, spike_train_driver, l_param, random_events); % create embeddings
TE_rate_1_2 = function_TE_rate(Cy, Jy, Cu_y, Ju_y, lambda_av_1_2, k_global);

%% MI rate
Cu_x = Cx;
Cu_y = Cy;
Ju_x = Jx; % or can be Jy

tau  = random_events(end, :) - random_events(1, :); % last spike, total recording duration!
lambda_av_u = Nu/tau; % average firing rate
MI_rate = function_MI_rate(Cu_x, Ju_x, Cu_y, lambda_av_u, k_global);

%% generate surrogates

ISI_1 = diff(spike_train1);  ISI_2 = diff(spike_train2);

isiSurr_1 = spiSeMe_surrogate_jodi(ISI_1, 'M', n_surr);
isiSurr_2 = spiSeMe_surrogate_jodi(ISI_2, 'M', n_surr);

TE_21_surr = []; TE_12_surr = []; MI_surr = [];

for k = 1 : size(isiSurr_1 ,2)
    %----------------------------------------------------------------------
    ISI_surr_1 =  isiSurr_1(:,k);
    spike_train_surr_1 = [0];
    for i = 2 : numel(ISI_surr_1)
        spike_train_surr_1(i) = ISI_surr_1(i) + spike_train_surr_1(i-1);
    end
    spike_train_surr_1 = spike_train_surr_1';
    %----------------------------------------------------------------------
    ISI_surr_2 =  isiSurr_2(:,k);
    spike_train_surr_2 = [0];
    for i = 2 : numel(ISI_surr_2)
        spike_train_surr_2(i) = ISI_surr_2(i) + spike_train_surr_2(i-1);
    end
    spike_train_surr_2 = spike_train_surr_2';
    
    %----------------------------------------------------------------------
    % TE rate 2 --> 1
    spike_train_target = spike_train_surr_1;
    spike_train_driver = spike_train_surr_2;
    
    [Cx, Jx, Cu_x, Ju_x] = function_embedding_vectors(spike_train_target, spike_train_driver, l_param, random_events); % create embeddings
    TE_21_surr_temp = function_TE_rate(Cx, Jx, Cu_x, Ju_x, lambda_av_2_1, k_global);
    TE_21_surr = [TE_21_surr; TE_21_surr_temp];
    
    %----------------------------------------------------------------------
    % TE rate 1 --> 2
    spike_train_target = spike_train_surr_2;
    spike_train_driver = spike_train_surr_1;
    
    [Cy, Jy, Cu_y, Ju_y] = function_embedding_vectors(spike_train_target, spike_train_driver, l_param, random_events); % create embeddings
    TE_12_surr_temp = function_TE_rate(Cy, Jy, Cu_y, Ju_y, lambda_av_1_2, k_global);
    TE_12_surr = [TE_12_surr; TE_12_surr_temp];
    
    %----------------------------------------------------------------------
    % MI rate 
    MI_surr_temp = function_MI_rate(Cu_x, Ju_x, Cu_y, lambda_av_u, k_global);
    MI_surr = [MI_surr; MI_surr_temp];
    
end % k, number of surrogates!


%% Visulization
figure(1);

row = 2;
col = 3;

subplot(row,col, [1:3]);
plot(spike_train1, 1*ones(1, numel(spike_train1)),'k.'); hold on; % train 1
plot(spike_train2, 2*ones(1, numel(spike_train2)), 'k.'); hold off; ylim([0 3]);
xlabel('Time [s]', 'FontName','Times New Roman','FontSize', 12);
yticks([1 2]); yticklabels({'Spike train 1','Spike train 2'}); ytickangle(45);
title('Coupled spike trains along direction 1-->2');
set(gca, 'FontName','Times New Roman','FontSize', 8);
axis tight;
ylim([0 3]);
%--------------------------------------------------------------------------
subplot(row,col, 4); % TE rate + surrogates distribution
plot(0.5, MI_rate, 'ro', 'MarkerFaceColor', 'r'); hold on;
plot(1, MI_surr, 'r*', 'MarkerFaceColor', 'r'); hold off;
xticks([0.5 1]); xticklabels({'Original value', 'Surrogates distr.'}); xtickangle(25); xlim([0 1.5]);
title('dMI rate [nats/s]', 'FontSize', 8);
set(gca, 'FontName','Times New Roman','FontSize', 8);

%--------------------------------------------------------------------------
subplot(row,col, 5); % TE rate + surrogates distribution
plot(0.5, TE_rate_1_2, 'bo', 'MarkerFaceColor', 'b'); hold on;
plot(1, TE_12_surr, 'b*', 'MarkerFaceColor', 'b'); hold off;

xticks([0.5 1]); xticklabels({'Original value', 'Surrogates distr.'}); xtickangle(25); xlim([0 1.5]);
title('TE rate 1-->2 [nats/s]', 'FontSize', 8);
set(gca, 'FontName','Times New Roman','FontSize', 8);

%--------------------------------------------------------------------------
subplot(row,col, 6); % TE rate + surrogates distribution
plot(0.5, TE_rate_2_1, 'go', 'MarkerFaceColor', 'g'); hold on;
plot(1, TE_21_surr, 'g*', 'MarkerFaceColor', 'g'); hold off;
xticks([0.5 1]); xticklabels({'Original value', 'Surrogates distr.'}); xtickangle(25); xlim([0 1.5]);
title('TE rate 2-->1[nats/s]', 'FontSize', 8);
set(gca, 'FontName','Times New Roman','FontSize', 8);

set(gcf, 'Color', 'w');
















