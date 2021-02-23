# TEMI


Here is demonstrated the framework that defines the mutual information rate (MIR) and the transfer entropy rate (TER) for two point processes X and Y, showing that the MIR between X and Y can be decomposed as the sum of the TER along the directions X → Y and Y → X, as proposed in the paper "An Information-Theoretic Framework to Measure the Dynamic Interaction between Neural Spike Trains" by G. Mijatovic, Y. Antonacci, T. Loncar Turukalo, L. Minati and L. Faes, 2021. 

The framework presents an information-theoretic approach for the model-free, continuous-time estimation of both undirected, symmetric (MIR) and directed, (Granger) causal (TER) interactions between spike trains (see REF0).

The TEMI toolbox contains functions for:

1. creating embeddings at target and joint events: function_embedding_vectors.m
2. estimation of transfer entropy rate: function_TE_rate.m 
4. assessment of measure' statistical significance based on surrogate data analysis: spiSeMe_surrogate_jodi.m (this function is part of the SpiSeMe package, see REF1)

and 

5. script demo.m to demonstrate the simulation of coupled spike trains denoted as spike train 1 (process X) and spike train 2 (process Y) where unidirectional interaction X → Y is simulated setting τ=δ= 0.1s (see the simulation protocol explained in Sect. III in the main article). In this simulation, both TER and MIR estimators are applied.

References:
REF0: Shorten, D., Spinney, R., & Lizier, J. (2020). Estimating Transfer Entropy in Continuous Time Between Neural Spike Trains or Other Event-Based Data. bioRxiv.
REF1: Ricci L, Castelluzzo M, Minati L, Perinelli A  (2019): "Generation  of  surro-gate event sequences via joint distribution of successive inter-event intervals", Chaos: An Interdisciplinary Journal of Nonlinear Science 29(12):121102


