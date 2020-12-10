function spikesPoisson = generatePoisson(rate, T)

% Generate spike train with Poisson distribution

% Input parameters: rate [number of spikes/s] and duration T [s]
% Output parameter: Poisson spike train

dt = 0.001; % % resolution of spike train in s, 1 ms

spikesPoisson = [];
for t = 0: dt : T
    if rate * dt >= rand(1)
        spikesPoisson = [spikesPoisson; t];
    end
end





