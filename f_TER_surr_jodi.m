function TER_surr_jodi = f_TER_surr_jodi(target, driver, num_surr, Nu, l_param, k_global)

ISIs_target = spiSeMe_surrogate_jodi(diff([0; target]), 'M', num_surr, 'verbose', false); % intervals
ISIs_driver = spiSeMe_surrogate_jodi(diff([0; driver]), 'M', num_surr, 'verbose', false); % intervals
spike_trains_s_t = cumsum(ISIs_target); % spikes
spike_trains_s_d = cumsum(ISIs_driver); % spikes

TER_surr_jodi = [];
for ns = 1 : num_surr
    
    spike_train_target = spike_trains_s_t(:,ns);
    spike_train_driver = spike_trains_s_d(:,ns);
    
    %% random events inside of each interval
    random_events = [];
    for st = 1 : numel(spike_train_target)-1
        start_lim = spike_train_target(st);
        end_lim = spike_train_target(st+1);
        random_events = [random_events; (end_lim-start_lim).*rand(Nu,1) + start_lim];
    end
    
    tau  = spike_train_target(end, :) - spike_train_target(1, :); % total recording duration!
    lambda = numel(spike_train_target)/tau; % average firing rate
    
    [Cx, Jx, Cu_x, Ju_x] = function_embedding_vectors(spike_train_target, spike_train_driver, l_param, random_events); % create embeddings
    TER_surr_temp = function_TE_rate(Cx, Jx, Cu_x, Ju_x, lambda, k_global);
    TER_surr_jodi = [TER_surr_jodi; TER_surr_temp];
end




