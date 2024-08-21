# Quantum anomalous Floquet phases
 
This code belongs to the manuscript "Quantum origin of anomalous Floquet phases in cavity-QED materials" [(arXiv:2312.10141)](https://arxiv.org/abs/2312.10141). 
It contains all the necessary simulations to reproduce the figures within the article. 

In fig1_energyspectrum.ipynb we compute and plot the energy spectrum of the total system, consisting on a SSH chain with photon-assisted hopping.

In fig2_topologicalinvariants.ipynb, we obtain the topological invariant for the total system, as well as the different contributions 
(i.e., from the single-particle system and the light-matter resonance). 

In fig3and4_dynamics.ipynb, we calculate the dynamics of the system for a given initial state, for different parameter choices. 
We also compute the Fourier transform of the result, to identify the main frequencies contributing to the complex dynamics. 

All the auxiliary functions are found in the scripts hamiltonians.py, dynamics.py, utils.py (for simulation and numerical results),
and line_gradient_colors.py (for plotting).