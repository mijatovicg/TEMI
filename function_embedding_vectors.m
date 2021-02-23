function [Cx, Jx, Cu, Ju] = function_embedding_vectors(spike_train_target, spike_train_driver, l_param, random_events)

% This function creates history and joint embeddings at target and random events

% Input parameters:
% two spike trains representing target and driver processes: spike_train_target and spike_train_driver
% embedding length: l_param
% stream of a randomly placed events: random events

% Output parameters:
% target history embeddings built at target events: Cx
% target history embeddings built at random events: Cu
% joint history embeddings built at target events: Jx
% joint history embeddings built at random events: Ju

%% input parsing & validation
if (~iscolumn(spike_train_target))
    spike_train_target = spike_train_target';
end
if (~iscolumn(spike_train_driver))
    spike_train_driver = spike_train_driver';
end

%% history embeddings of target and driver spike trains in respect to target events
Cx = []; Jx = [];
for i = l_param+1 : numel(spike_train_target) % start in the way to have past of target process for sure
    
    temp_spike = spike_train_target(i, :); % target event
    
    if numel(find(spike_train_driver < temp_spike))>= l_param  %  only if past of driver exists
        
        % spike_train_target
        target_past = spike_train_target(find(spike_train_target < temp_spike)); % past of target, strictly less than target event(temp_spike)
        temp_target = [];
        temp_target = [temp_target; temp_spike; target_past(end:-1:end-l_param+1, :)]; % take past of target --> l_param+1 spikes
        temp_target = temp_target(end:-1:1, :); % take past of target --> l_param+1 spikes, put in desired order
        
        % spike_train_driver
        driver_past = spike_train_driver(find(spike_train_driver < temp_spike)); % past of driver, strictly less than target event(temp_spike)
        temp_driver = [];
        temp_driver = [temp_driver; temp_spike; driver_past(end:-1:end-l_param+1, :)]; % take past of driver --> l_param+1 spikes
        temp_driver = temp_driver(end:-1:1, :); % take past of driver --> l_param+1 spikes, put in desired order
        
        % matrices of embedding vectors
        Cx = [Cx; diff(temp_target')];
        Jx = [Jx; diff(temp_target') diff(temp_driver')];
    end
end

%% history embeddings of target and driver spike trains in respect to random events
Cu = []; Ju = [];
for i = 1 : numel(random_events)
    
    temp_spike = random_events(i, :); % random event
    
    target_past = spike_train_target(find(spike_train_target < temp_spike)); % past of target, stricly less than random event
    driver_past = spike_train_driver(find(spike_train_driver < temp_spike)); % past od driver, stricly less than random event
    
    if (numel(target_past)>= l_param) && (numel(driver_past)>= l_param) % if both pasts exist
        
        % spike_train_target
        temp_target = [];
        temp_target = [temp_target; temp_spike; target_past(end:-1:end-l_param+1, :)]; % take past of target --> l_param+1 spikes
        temp_target = temp_target(end:-1:1, :); % take past of target --> l_param+1 spikes, put in desired order
        
        % spike_train_driver
        temp_driver = [];
        temp_driver = [temp_driver; temp_spike; driver_past(end:-1:end-l_param+1, :)]; % take past of driver --> l_param+1 spikes
        temp_driver = temp_driver(end:-1:1, :); % take past of driver --> l_param+1 spikes, put in desired order
        
        % matices of embedding vectors
        Cu = [Cu; diff(temp_target')];
        Ju = [Ju; diff(temp_target') diff(temp_driver')];
        
    end
end
