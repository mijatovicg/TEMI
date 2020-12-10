# TEMI


Here is demonstrated the framework for the transfer entropy (TE) rate and  dynamical mutual information (dMI) estimation, proposed in the paper 
"An Information−Theoretic Framework to Measure the Dynamic Interaction between Neural Spike Trains" by G. Mijatovic, Y. Antonacci, T. Loncar Turukalo, L. Minati and L. Faes,2020. The framework presents an information-theoretic approach for the model-free, continuous-time estimation of of both undirected, symmetric (dMI) and directed, causal (TE) interactions between spike trains.

The TEMI toolbox contains functions for:

1. creating embeddings at target and joint events: function_embedding_vectors.m
2. estimation of transfer entropy rate: function_TE_rate
3. estimation of dynamic mutual information rate: function_MI_rate
4. assessment of measure' statistical significance based on surrogate data analysis: spiSeMe_surrogate_jodi.m (this function is part of the SpiSeMe package see REF0)

and 

5. script demo.m to demonstrate the simulation of coupled spike trains denoted as spike train 1 (process X) and 2 (process Y) where unidirectional interaction X→Y is simulated setting τ=δ= 0.1s (see the simulation protocol explained in Sect. III in the main article).

REF0: Ricci L, Castelluzzo M, Minati L, Perinelli A  (2019): "Generation  of  surro-gate event sequences via joint distribution of successive inter-event intervals", Chaos: An Interdisciplinary Journal of Nonlinear Science 29(12):121102
